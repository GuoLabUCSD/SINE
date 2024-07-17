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
args_path.add_argument('-j', '--junctions_table', type=str, required=True, help='Full junctions table generated from get_wt.sh')
args_path.add_argument('-a', '--affinity_scores', type=str, required=True, help='Directory with allele affinity scores. Ie) Output files from xls_parser.sh')
args_path.add_argument('-m', '--mapping_directory', type=str, required=True, help='Directory containing patient files with the mapping ID number to netmhcpan')
args_path.add_argument('-o', '--output_directory', type=str, required=True, help='Location to store phbr score results and plots')
args_path.add_argument('-t', '--paired_table', type=str, required=True, help='Location of the ASE junctions table with their corresponding paired WT junctions')

args = parser.parse_args()

def aggregate_ase(PHBR_col):
    final_PHBR_column_to_add = []
    for event_scores in PHBR_col:
        phbr_average = (sum(event_scores) / len(event_scores))
        final_PHBR_column_to_add.append(phbr_average)
    return final_PHBR_column_to_add

def aggregate_wt(pep_col, PHBR_col):
    final_peptide_column_to_add = []
    final_PHBR_column_to_add = []
    for event, event_scores in zip(pep_col, PHBR_col):
        if event == 'nan':
            final_PHBR_column_to_add.append('nan')
        else:
            phbr_average = (sum(event_scores) / len(event_scores))
            final_PHBR_column_to_add.append(phbr_average)
    return final_PHBR_column_to_add

# Load junctions we generated neopeptides for
junc_df = pd.read_csv('{}'.format(args.junctions_table), sep='\t')
samples_only_df = junc_df.drop_duplicates(subset = 'patient')
tums = samples_only_df[samples_only_df['junction_category'] == 'ASE_junc_of_interest']
tum_patients = [i for i in tums['patient']]
norms = samples_only_df[samples_only_df['junction_category'] == 'other']

# Gather NETMHCpan parsed output affinity files for each patient
total_df = pd.DataFrame()
for patient in tum_patients:
    aff_path = '{}/{}.output'.format(args.affinity_scores, patient)
    maf_path = '{}/{}.tsv'.format(args.mapping_directory, patient)

    try:
        # load file to map NetMHCpan short IDs to the actual junction names and isoform
        pat_df = pd.read_csv(maf_path, sep='\t', usecols=['junction','isoform','NetMHCpan_short_ID', 'junc_aa_pos'])

        # load patient affinity file
        pat_aff = pd.read_csv(aff_path, sep='\t', index_col=0)

        # join
        pat_aff = pat_aff.join(pat_df).set_index('NetMHCpan_short_ID')

        #Drop Peptide Columns to assign back later
        #print(pat_aff)
        peptide_df = pat_aff[['junction','isoform', 'peptide', 'junc_aa_pos']]
        #peptide_df = pat_aff.drop(['H-2-Kb', 'H-2-Db'], axis=1)
        #print(peptide_df)
        pat_aff = pat_aff.drop(['peptide', 'isoform', 'junc_aa_pos', 'junction'], axis=1)

        # take harmonic mean as in Marty et al. 2017
        # find the highest binding allele for each event
        pat_aff['Best_Allele'] = pat_aff.T.idxmin()
        allele_df = pat_aff[['Best_Allele']]
        pat_aff = pat_aff.drop(['Best_Allele'], axis = 1)
        pat_aff = pat_aff.apply(harmonic_mean, axis=1)
        pat_aff = pd.DataFrame(pat_aff, columns=['PHBR_scores']).reset_index()

        #Join Peptides back to the df
        peptide_df = peptide_df.reset_index()
        #print(peptide_df)
        pat_aff = pd.merge(pat_aff, peptide_df, on = 'NetMHCpan_short_ID', how = 'left')
        pat_aff = pd.merge(pat_aff, allele_df, on = 'NetMHCpan_short_ID', how = 'left')

        clean_pat_aff = pd.DataFrame(columns = ['NetMHCpan_short_ID', 'PHBR_scores', 'peptide', 'isoform', 'junc_aa_pos', 'junction', 'Best_Allele'])
        new_index = 0
        for index, row in pat_aff.iterrows():
            clean_pos = literal_eval(row['junc_aa_pos'])
            result = re.finditer(r'(?=({}))'.format(row['peptide']), row['isoform'])
            for match in result:
                if clean_pos[0] in range(match.start(1), match.end(1)) and clean_pos[1] in range(match.start(1), match.end(1)):
                    clean_pat_aff.loc[new_index] = row
                    new_index = new_index + 1
                    continue
        
        #Merge with patient specific values from original table
        #Make wild-type rows show the patient from which HLA-types were checked against
        pat_junc_df = junc_df[(junc_df['patient'] == patient) | (junc_df['junction_category'] == 'other')]
        pat_junc_df = pat_junc_df.assign(patient = patient)

        # add manual annotations
        clean_pat_aff = pd.merge(clean_pat_aff, pat_junc_df[['junction_category','junction','chr_intron_start','chr_intron_end','num_uniq_reads_at_junc','patient']],
                           on='junction', how='left')
        clean_pat_aff = clean_pat_aff.drop(['NetMHCpan_short_ID'], axis = 1)

        clean_pat_aff['PHBR_scores'] = pd.to_numeric(clean_pat_aff['PHBR_scores'])


        clean_pat_aff = clean_pat_aff[clean_pat_aff['PHBR_scores'] == clean_pat_aff.groupby('junction')['PHBR_scores'].transform(min)]


        clean_pat_aff = clean_pat_aff.groupby(['PHBR_scores', 'peptide', 'junction', 'Best_Allele', 'junction_category', 'chr_intron_start', 'chr_intron_end', 
                                               'num_uniq_reads_at_junc', 'patient'])[['isoform', 'junc_aa_pos']].agg(list).reset_index()

        total_df = total_df.append(clean_pat_aff)

    except:
        continue

total_df = total_df.drop(columns = ['chr_intron_start', 'chr_intron_end'], axis = 1)
wt_temp = total_df.loc[total_df['junction_category'] == 'other']
total_df = total_df[total_df['junction_category'] != 'other']
wt_temp = wt_temp.drop_duplicates(subset = ['PHBR_scores', 'peptide', 'junction', 'patient'])
total_df = total_df.append(wt_temp)

#Split the combined dataframe of ASE and WT into separate frames
ase_frame = total_df[total_df['junction_category'] == 'ASE_junc_of_interest']
ase_frame = ase_frame.drop(columns = ['junction_category'], axis = 1)
wt_frame = total_df[total_df['junction_category'] == 'other']
wt_frame = wt_frame.drop(columns = ['junction_category'], axis = 1)

all_ases = ase_frame.groupby('junction')
all_wts = wt_frame.groupby('junction')

#Group both dataframes so there is just one junction per row and each cell contains a list of all values found in each patient
final_ase_df = pd.DataFrame(columns = ['ASE_PHBR_Score', 'ASE_peptide', 'ASE_junction', 'ASE_#Reads', 'Tumor_samples_with_junction', '%samples_with_ASE_junction', 'ASE_isoforms', 'ASE_junction_index', 'ASE_best_allele'])
final_wt_df = pd.DataFrame(columns = ['WT_PHBR_Score', 'WT_peptide', 'WT_junction', 'WT_#Reads', 'Tumor_Sample_Map_to_WT_PHBR', 'WT_isoforms', 'WT_junction_index', 'WT_best_allele'])

row_index = 0
for name, group in all_ases:
    scores = group['PHBR_scores'].tolist()
    peptides = group['peptide'].tolist()
    alleles = group['Best_Allele'].to_list()
    junctions = name
    reads = group['num_uniq_reads_at_junc'].to_list()
    samples = group['patient'].to_list()
    all_isos = group['isoform'].to_list()
    perc = (len(samples) / len(tums)) * 100
    all_positions = group['junc_aa_pos'].to_list()
    final_ase_df.loc[row_index] = [scores, peptides, junctions, reads, samples, perc, all_isos, all_positions, alleles]
    row_index = row_index + 1

row_index = 0
for name, group in all_wts:
    scores = group['PHBR_scores'].tolist()
    peptides = group['peptide'].tolist()
    alleles = group['Best_Allele'].to_list()
    junctions = name
    reads = group['num_uniq_reads_at_junc'].to_list()
    samples = group['patient'].to_list()
    all_isos = group['isoform'].to_list()
    all_positions = group['junc_aa_pos'].to_list()
    final_wt_df.loc[row_index] = [scores, peptides, junctions, reads, samples, all_isos, all_positions, alleles]
    row_index = row_index + 1

final_ase_df = final_ase_df.set_index('ASE_junction', drop = False)
final_wt_df = final_wt_df.set_index('WT_junction', drop = False)

#Load Pairing Table
pairing_table = pd.read_csv('{}'.format(args.paired_table), sep = '\t')
pairing_table['wt'] = pairing_table['wt'].fillna(0)

#For every row in the pairing table, merge the corresponding ase row with the corresponding wt row
my_final_df = pd.DataFrame(columns = ['ASE_junction', 'ASE_peptide', 'ASE_PHBR_Score', 'ASE_#Reads', 'Tumor_samples_with_junction', '%samples_with_ASE_junction', 'ASE_isoforms', 'ASE_junction_indices', 'ASE_best_alleles', 'WT_junction', 'WT_peptide', 'WT_PHBR_Score', 'WT_#Reads', 'Tumor_Sample_Map_to_WT_PHBR', 'WT_isoforms', 'WT_junction_indices', 'WT_best_alleles'])

index_val = 0
for index, row in pairing_table.iterrows():
    if row['ase'] in final_ase_df.index:
        if row['wt'] == 0 or row['wt'] not in final_wt_df.index:
            x = final_ase_df.loc[row['ase']]
            val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12, val13, val14, val15, val16, val17 = x['ASE_junction'], x['ASE_peptide'], x['ASE_PHBR_Score'], x['ASE_#Reads'], x['Tumor_samples_with_junction'], x['%samples_with_ASE_junction'], x['ASE_isoforms'], x['ASE_junction_index'], x['ASE_best_allele'], 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan'
            my_final_df.loc[index_val] = val1,val2,val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,val15,val16,val17
            index_val = index_val + 1
        else:
            x = final_ase_df.loc[row['ase']]
            y = final_wt_df.loc[row['wt']]
            val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12, val13, val14, val15, val16, val17 = x['ASE_junction'], x['ASE_peptide'], x['ASE_PHBR_Score'], x['ASE_#Reads'], x['Tumor_samples_with_junction'], x['%samples_with_ASE_junction'], x['ASE_isoforms'], x['ASE_junction_index'], x['ASE_best_allele'],  y['WT_junction'], y['WT_peptide'], y['WT_PHBR_Score'], y['WT_#Reads'], y['Tumor_Sample_Map_to_WT_PHBR'], y['WT_isoforms'], y['WT_junction_index'], y['WT_best_allele']
            my_final_df.loc[index_val] = val1,val2,val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,val15,val16,val17
            index_val = index_val + 1
    else:
        continue
my_final_df = my_final_df.drop(columns = ['WT_#Reads', 'WT_isoforms', 'WT_junction_indices', 'ASE_junction_indices'])
raw_df = my_final_df.copy()

ASE_PHBR_col = my_final_df['ASE_PHBR_Score'].to_list()
WT_pep_col = my_final_df['WT_peptide'].to_list()
WT_PHBR_col = my_final_df['WT_PHBR_Score'].to_list()
b = my_final_df['ASE_peptide'].to_list()
t = my_final_df['ASE_best_alleles'].to_list()
c = my_final_df['ASE_PHBR_Score'].to_list()
best_pairs_2frame = []
for i,j,k in zip(b,t,c):
    min_position = k.index(min(k))
    best_pairs_2frame.append([i[min_position], j[min_position]])
b = my_final_df['WT_peptide'].to_list()
t = my_final_df['WT_best_alleles'].to_list()
c = my_final_df['WT_PHBR_Score'].to_list()
best_pairs_2frameWT = []
for i,j,k in zip(b,t,c):
    min_position = k.index(min(k))
    best_pairs_2frameWT.append([i[min_position], j[min_position]])
x = aggregate_ase(ASE_PHBR_col)
y = aggregate_wt(WT_pep_col, WT_PHBR_col)
my_final_df['ASE_PHBR_Score'] = x
my_final_df['WT_PHBR_Score'] = y
my_final_df = my_final_df.drop(columns = ['ASE_#Reads', 'ASE_isoforms', 'Tumor_Sample_Map_to_WT_PHBR'])
my_final_df['Best_ASE_Peptide_HLA_Pair'] = best_pairs_2frame
my_final_df['Best_WT_Peptide_HLA_Pair'] = best_pairs_2frameWT
my_final_df = my_final_df[['ASE_junction', 'Best_ASE_Peptide_HLA_Pair', 'ASE_peptide', 'ASE_PHBR_Score', 'Tumor_samples_with_junction', '%samples_with_ASE_junction', 'ASE_best_alleles', 'WT_junction', 'WT_peptide', 'WT_PHBR_Score', 'WT_best_alleles']]
my_final_df = my_final_df.rename(columns = {'ASE_PHBR_Score': 'Average_ASE_PHBR_Score(s)', 'WT_PHBR_Score': 'Average_WT_PHBR_Score(s)'}) 

#Loop through all junction to gene mapping files in each sample directory (./isoforms/tumors)
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
                if row['Junction'] not in all_juncGenes:
                    all_juncGenes[row['Junction']] = row['Symbol']
                else:
                    continue

raw_df['Symbol'] = raw_df['ASE_junction'].map(all_juncGenes)
my_final_df['Symbol'] = my_final_df['ASE_junction'].map(all_juncGenes)

#Remove rows with a junction/gene strand mismatch
#First loop through all of the strand info files and get a dictionary of junction / gene mismatches
#Drop rows where both the junction and the gene are in that dictionary
path = '{}'.format(args.isoforms_dir)
mismatch_dict = {}

for sample_type in os.listdir(path):
    if sample_type == 'tumors':
        path = os.path.join(path, 'tumors')
        for sample_dir in os.listdir(path):
            strand_file_path = os.path.join(path, sample_dir, 'intermediate_results', 'strand_flags.txt')
            strand_file = pd.read_csv(strand_file_path, sep = '\t')
            strand_file = strand_file.loc[strand_file['Flag'] == 'MISMATCH']
            for index, row in strand_file.iterrows():
                if row['Junction'] not in mismatch_dict:
                    mismatch_dict[row['Junction']] = row['Symbol']

output_df = my_final_df.copy()
for index, row in my_final_df.iterrows():
    pair = (row['ASE_junction'], row['Symbol'])
    for key, value in mismatch_dict.items():
        if (key, value) == pair:
            output_df = output_df[(output_df['ASE_junction'] != key) & (output_df['Symbol'] != value)]

for index, row in output_df.iterrows():
    b_pep = row['Best_ASE_Peptide_HLA_Pair'][0]
    b_hla = row['Best_ASE_Peptide_HLA_Pair'][1]
    samples = row['Tumor_samples_with_junction']
    for peptide, allele, sample in zip(row['ASE_peptide'], row['ASE_best_alleles'], row['Tumor_samples_with_junction']):
        if peptide == b_pep and allele == b_hla:
            file = pd.read_csv('{}/{}.output'.format(args.affinity_scores, sample), sep = '\t')
            file = file.loc[file['peptide'] == peptide].reset_index(drop = True)
            first_row = file.iloc[[0]]
            first_row = first_row.drop(columns = ['Unnamed: 0', 'peptide'])
            value = first_row.min(axis = 1)
            row['Best_ASE_Peptide_HLA_Pair'].append(float(value))
            break

output_df.to_csv('{}/fin_results.txt'.format(args.output_directory), sep='\t', index=False)
raw_df.to_csv('{}/raw_results.txt'.format(args.output_directory), sep='\t', index=False)

#Make Plots
#def plot_boxplots(df, x, y=None, y_list=None, order=None, ylim_dict=None, figsize=(7,4), hue=None, palette='Set1',
#        savepath=None, suptitle=None, comparisons_correction=None, plot_significant_only=True):
#        if y and y_list:
#                print('ERROR. Cannot have y and y_list at the same time')
#                return
#
#
#        plt.figure(figsize=figsize)
#        if order is None:
#                order = df[x].unique()
#
#
#        if y_list:
#
#                i=1
#
#                for y in y_list:
#                        plt.subplot(1,len(y_list),i)
#                        ax = sns.boxplot(x=x, y=y, data=df, order=order, palette=palette)
#
#                        if ylim_dict:
#                                if y in ylim_dict:
#                                        plt.ylim(ylim_dict[y])
#
#                        add_stat_annotation(ax=ax, x=x, y=y, data=df, test='Mann-Whitney', text_format='simple', order=order, plot_significant_only=plot_significant_only,
#                                                                box_pairs=create_pairs(order), verbose=0, loc='outside', comparisons_correction=comparisons_correction)
#                        plot_medians_in_boxplot(ax=ax, df=df, x=x, y=y)
#                        counts = df[x].value_counts()
#                        ax.set_xticklabels(['{}\n({})'.format(x.get_text(),counts.loc[x.get_text()]) for x in ax.get_xticklabels()])
#                        plt.xlabel('')
#                        i+=1
#                plt.tight_layout()
#
#                if suptitle:
#                        plt.suptitle(suptitle, y=1.05)
#
#        else:
#                if hue:
#                        hue_order = df[hue].unique()
#
#                        ax = sns.boxplot(x=x, y=y, data=df, order=order, hue=hue, hue_order=hue_order, palette=palette)
#
#                        if ylim_dict:
#                                if y in ylim_dict:
#                                        plt.ylim(ylim_dict[y])
#
#                        add_stat_annotation(ax=ax, x=x, y=y, data=df, test='Mann-Whitney', text_format='simple', order=order, hue=hue, hue_order=hue_order,
#                                                                box_pairs=create_hue_pairs(order, hue_order), verbose=0, loc='outside', comparisons_correction=comparisons_correction)
#                        if len(hue_order)>2:
#                                print('Can only handle plotting medians in hue pairs of 2.')
#                        else:
#                                plot_medians_in_boxplot(ax=ax, df=df, x=x, y=y, hue_left_right=hue_order, hue=hue)
#                                counts = df[[x,hue]].value_counts()
#                                ax.set_xticklabels([f'{x.get_text()}\n({counts.loc[x.get_text()][hue_order[0]]},{counts.loc[x.get_text()][hue_order[1]]})' for x in ax.get_xticklabels()])
#                        plt.xlabel('')
#                else:
#                        ax = sns.boxplot(x=x, y=y, data=df, order=order, palette=palette)
#
#                        if ylim_dict:
#                                if y in ylim_dict:
#                                        plt.ylim(ylim_dict[y])
#
#                        add_stat_annotation(ax=ax, x=x, y=y, data=df, test='Mann-Whitney', text_format='simple', order=order,
#                                                                box_pairs=create_pairs(order), verbose=0, loc='outside', comparisons_correction=comparisons_correction)
#                        plot_medians_in_boxplot(ax=ax, df=df, x=x, y=y)
#                        counts = df[x].value_counts()
#                        ax.set_xticklabels(['{}\n({})'.format(x.get_text(),counts.loc[x.get_text()]) for x in ax.get_xticklabels()])
#                        plt.xlabel('')
#        if savepath:
#                print(savepath)
#                plt.savefig(savepath, bbox_inches='tight')
#                plt.show()
#        return
#
#plot_boxplots(x='junction_category', y='PHBR_scores', df=total_df.dropna(), figsize=(3,3), ylim_dict={'PHBR_scores':(-1,21)}, 
#              order=['other','ASE_junc_of_interest'], savepath='{}/ase_vs_other.png'.format(args.output_directory))
