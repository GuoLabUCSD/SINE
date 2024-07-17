import argparse
import pandas as pd
import numpy as np
import os
from scipy.stats import hmean as harmonic_mean
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from ast import literal_eval
import re
import statistics
#import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
#import seaborn as sns
#sns.set_style('white')
#sns.set_style('ticks')
#from statannot import add_stat_annotation
#import matplotlib as mpl
#mpl.rcParams['pdf.fonttype'] = 42
#from plot_medians_in_boxplot import plot_medians_in_boxplot
#from stat_annot_utils import create_pairs, create_hue_pairs

parser = argparse.ArgumentParser()
args_path = parser.add_argument_group("Inputs:")

args_path.add_argument('-i', '--isoforms_dir', type=str, required=True, help='Directory of isoforms generated from main.sh')
args_path.add_argument('-a', '--affinity_scores', type=str, required=True, help='Directory with allele affinity scores. Ie) Output files from xls_parser.sh')
args_path.add_argument('-m', '--mapping_directory', type=str, required=True, help='Directory containing patient files with the mapping ID number to netmhcpan')
args_path.add_argument('-o', '--output_directory', type=str, required=True, help='Location to store phbr score results and plots')

args = parser.parse_args()

def aggregate_ase(PHBR_col):
    final_PHBR_column_to_add = []
    for event_scores in PHBR_col:
        phbr_average = (sum(event_scores) / len(event_scores))
        final_PHBR_column_to_add.append(phbr_average)
    return final_PHBR_column_to_add

results_path = '{}'.format(args.affinity_scores)
tum_patients = []
for file in os.listdir(results_path):
    file = file.split('.')
    tum_patients.append(file[0])

total_df = pd.DataFrame()
for patient in tum_patients:
    aff_path = '{}/{}.output'.format(args.affinity_scores, patient)
    maf_path = '{}/{}.tsv'.format(args.mapping_directory, patient)
    
    # load file to map NetMHCpan short IDs to the actual junction names and isoform
    pat_df = pd.read_csv(maf_path, sep='\t', usecols=['insertion','isoform','NetMHCpan_short_ID', 'ins_aa_pos'])

    # load patient affinity file
    pat_aff = pd.read_csv(aff_path, sep='\t', index_col=0)

    # join
    pat_aff = pat_aff.join(pat_df).set_index('NetMHCpan_short_ID')

    #Drop Peptide Columns to assign back later
    #print(pat_aff)
    peptide_df = pat_aff[['insertion','isoform', 'peptide', 'ins_aa_pos']]
    #peptide_df = pat_aff.drop(['H-2-Kb', 'H-2-Db'], axis=1)
    #print(peptide_df)
    pat_aff = pat_aff.drop(['peptide', 'isoform', 'ins_aa_pos', 'insertion'], axis=1)

    # take harmonic mean as in Marty et al. 2017
    # find the highest binding allele for each event
    pat_aff['Best_Allele'] = pat_aff.T.idxmin()
    allele_df = pat_aff[['Best_Allele']]
    pat_aff = pat_aff.drop(['Best_Allele'], axis = 1)
    pat_aff = pat_aff.apply(harmonic_mean, axis=1)
    pat_aff = pd.DataFrame(pat_aff, columns=['PHBR_scores']).reset_index()

    #Join Peptides and best alleles back to the df
    peptide_df = peptide_df.reset_index()
    #print(peptide_df)
    pat_aff = pd.merge(pat_aff, peptide_df, on = 'NetMHCpan_short_ID', how = 'left')
    pat_aff = pd.merge(pat_aff, allele_df, on = 'NetMHCpan_short_ID', how = 'left')

    clean_pat_aff = pd.DataFrame(columns = ['NetMHCpan_short_ID', 'PHBR_scores', 'peptide', 'isoform', 'ins_aa_pos', 'insertion', 'Best_Allele'])
    new_index = 0
    for index, row in pat_aff.iterrows():
        #clean_pos = literal_eval(row['ins_aa_pos'])
        #result = re.finditer(r'(?=({}))'.format(row['peptide']), row['isoform'])
        #for match in result:
            #if clean_pos[0] in range(match.start(1), match.end(1)) and clean_pos[1] in range(match.start(1), match.end(1)):
        clean_pat_aff.loc[new_index] = row
        new_index = new_index + 1
                #continue

    # add manual annotations
    # if nan, insertion is not located on the spot of the isoform
    clean_pat_aff = clean_pat_aff.drop(['NetMHCpan_short_ID'], axis = 1)
    
    clean_pat_aff['patient'] = patient

    clean_pat_aff['PHBR_scores'] = pd.to_numeric(clean_pat_aff['PHBR_scores'])


    clean_pat_aff = clean_pat_aff[clean_pat_aff['PHBR_scores'] == clean_pat_aff.groupby('insertion')['PHBR_scores'].transform(min)]


    clean_pat_aff = clean_pat_aff.groupby(['PHBR_scores', 'peptide', 'insertion', 'Best_Allele', 
                                           'patient'])[['isoform', 'ins_aa_pos']].agg(list).reset_index()

    total_df = total_df.append(clean_pat_aff)

all_ases = total_df.groupby('insertion')

final_ase_df = pd.DataFrame(columns = ['ASE_PHBR_Score', 'ASE_peptide', 'ASE_insertion', 'Tumor_samples_with_insertion', '%samples_with_ASE_insertion', 'ASE_isoforms', 'ASE_insertion_index', 'ASE_best_allele'])

row_index = 0
for name, group in all_ases:
    scores = group['PHBR_scores'].tolist()
    peptides = group['peptide'].tolist()
    alleles = group['Best_Allele'].to_list()
    junctions = name
    samples = group['patient'].to_list()
    all_isos = group['isoform'].to_list()
    perc = (len(samples) / len(tum_patients)) * 100
    all_positions = group['ins_aa_pos'].to_list()
    final_ase_df.loc[row_index] = [scores, peptides, junctions, samples, perc, all_isos, all_positions, alleles]
    row_index = row_index + 1

final_ase_df = final_ase_df.set_index('ASE_insertion', drop = False)

my_final_df = final_ase_df.reset_index(drop = True)
my_final_df = my_final_df.rename(columns = {'ASE_best_allele':'ASE_best_alleles'})
raw_df = my_final_df.copy()

ASE_PHBR_col = my_final_df['ASE_PHBR_Score'].to_list()
b = my_final_df['ASE_peptide'].to_list()
t = my_final_df['ASE_best_alleles'].to_list()
c = my_final_df['ASE_PHBR_Score'].to_list()
best_pairs_2frame = []
for i,j,k in zip(b,t,c):
    min_position = k.index(min(k))
    best_pairs_2frame.append([i[min_position], j[min_position]])
x = aggregate_ase(ASE_PHBR_col)
my_final_df['ASE_PHBR_Score'] = x
my_final_df = my_final_df.drop(columns = ['ASE_isoforms'])
my_final_df['Best_ASE_Peptide_HLA_Pair'] = best_pairs_2frame
my_final_df = my_final_df[['ASE_insertion', 'Best_ASE_Peptide_HLA_Pair', 'ASE_peptide', 'ASE_PHBR_Score', 'Tumor_samples_with_insertion', '%samples_with_ASE_insertion', 'ASE_best_alleles']]
my_final_df = my_final_df.rename(columns = {'ASE_PHBR_Score': 'Average_ASE_PHBR_Score(s)'}) 

path = '{}'.format(args.isoforms_dir)
all_juncGenes = {}

for sample_type in os.listdir(path):
    if sample_type == 'tumors':
        path = os.path.join(path, 'tumors')
        for sample_dir in os.listdir(path):
            sample_genes = {}
            gene_file_path = os.path.join(path, sample_dir, 'intermediate_results', 'gene_mappings.txt')
            gene_file = pd.read_csv(gene_file_path, sep = '\t')
            for index, row in gene_file.iterrows():
                if row['ROI'] not in all_juncGenes:
                    all_juncGenes[row['ROI']] = row['Symbol']
                else:
                    continue

raw_df['Symbol'] = raw_df['ASE_insertion'].map(all_juncGenes)
my_final_df['Symbol'] = my_final_df['ASE_insertion'].map(all_juncGenes)

path = '{}'.format(args.isoforms_dir)
mismatch_dict = {}

for sample_type in os.listdir(path):
    if sample_type == 'tumors':
        path = os.path.join(path, 'tumors')
        for sample_dir in os.listdir(path):
            try:
                strand_file_path = os.path.join(path, sample_dir, 'intermediate_results', 'strand_flags.txt')
                strand_file = pd.read_csv(strand_file_path, sep = '\t')
                strand_file = strand_file.loc[strand_file['Flag'] == 'MISMATCH']
                for index, row in strand_file.iterrows():
                    if row['ROI'] not in mismatch_dict:
                        mismatch_dict[row['ROI']] = row['Symbol']
            except:
                continue

output_df = my_final_df.copy()
for index, row in my_final_df.iterrows():
    pair = (row['ASE_insertion'], row['Symbol'])
    for key, value in mismatch_dict.items():
        if (key, value) == pair:
            output_df = output_df[(output_df['ASE_junction'] != key) & (output_df['Symbol'] != value)]

for index, row in output_df.iterrows():
    b_pep = row['Best_ASE_Peptide_HLA_Pair'][0]
    b_hla = row['Best_ASE_Peptide_HLA_Pair'][1]
    samples = row['Tumor_samples_with_junction']
    for peptide, allele, sample in zip(row['ASE_peptide'], row['ASE_best_alleles'], row['Tumor_samples_with_insertion']):
        if peptide == b_pep and allele == b_hla:
            file = pd.read_csv('{}/{}.output'.format(args.affinity_scores, sample), sep = '\t')
            file = file.loc[file['peptide'] == peptide].reset_index(drop = True)
            first_row = file.iloc[[0]]
            first_row = first_row.drop(columns = ['Unnamed: 0', 'peptide'])
            value = first_row.min(axis = 1)
            row['Best_ASE_Peptide_HLA_Pair'].append(float(value))
            break

output_df = output_df.rename(columns = {'ASE_insertion': 'ASE_ROI', 'Tumor_samples_with_insertion': 'Tumor_samples_with_ROI', '%samples_with_ASE_insertion': '%samples_with_ASE_ROI'}) 
raw_df = raw_df.rename(columns = {'ASE_insertion': 'ASE_ROI', 'Tumor_samples_with_insertion': 'Tumor_samples_with_ROI', '%samples_with_ASE_insertion': '%samples_with_ASE_ROI'}) 
output_df.to_csv('{}/fin_results.txt'.format(args.output_directory), sep='\t', index=False)
raw_df.to_csv('{}/raw_results.txt'.format(args.output_directory), sep='\t', index=False)

