import argparse
import pandas as pd
import numpy as np
import os
from ast import literal_eval
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

parser = argparse.ArgumentParser()
args_path = parser.add_argument_group("Inputs:")

args_path.add_argument('-i', '--isoforms_dir', type=str, required=True, help='Directory of isoforms generated from main.sh')
args_path.add_argument('-f', '--output_directory_fastas', type=str, required=True, help='Location to store output fasta files')
args_path.add_argument('-m', '--output_directory_mapping', type=str, required=True, help='Location to store files mapping patient junctions to a netmhcpan ID')
args_path.add_argument('-o', '--netmhcpan_output_directory', type=str, required=True, help='Location to store netmhcpan output')
args_path.add_argument('-a', '--allele_table', type=str, required=True, help='Filepath to a tab-separated text file of MHC alleles to check peptide binding against for each sample')
args_path.add_argument('-s', '--patient_sample', type=str, required=True, help='Patient name to run netmhcpan on. Ie) Sample Directory name in the isoforms directory')
args_path.add_argument('-p', '--pipeline_directory', type=str, required=True, help='Directory to the pipeline')
args_path.add_argument('-n', '--net_two', type=str, required=False, help='Use NetMHCIIpan instead')

args = parser.parse_args()

# Helper Function

def create_peptides(mhc, short_seq, junc_aa_pos):
    '''Creates peptides spanning junction from buffered sequence'''
    kmer_list = [8,9,10,11] if mhc == 'I' else [15]
    peptides = []

    for kmer in kmer_list:
        for i in range(len(short_seq) - (kmer-1)):
            start = i
            end = i + kmer 
            # check if spans junctions
            if start <= junc_aa_pos[0] and junc_aa_pos[1] < end:
                peptides.append(short_seq[start:end])
    return peptides

#Concatenate the WT junction peptides found in the normals all together
path = '{}'.format(args.isoforms_dir)
norm_path = os.path.join(path, 'normals')
full_wt_peptides = pd.DataFrame(columns = ['isoform', 'isoform_num', 'junc_aa_pos', 'peptides', 'junction', 'frame_selection'])
for normal in os.listdir(norm_path):
    wt_peptides_path = os.path.join(norm_path, normal, 'intermediate_results', 'neopeptides.tsv')
    peps = pd.read_csv(wt_peptides_path, sep = '\t')
    full_wt_peptides = full_wt_peptides.append(peps)
full_wt_peptides = full_wt_peptides.drop_duplicates(subset = ['isoform', 'junc_aa_pos', 'peptides', 'junction', 'frame_selection'])
full_wt_peptides = full_wt_peptides.reset_index(drop = True)

#Concatenate this new master list of normal peptides to the neopeptides.tsv files (might as well save it as a new file) for each sample
tum_path = os.path.join(path, 'tumors')
for tumor in os.listdir(tum_path):
    tum_peptides_path = os.path.join(tum_path, tumor, 'intermediate_results', 'neopeptides.tsv')
    intermediate_dir = os.path.join(tum_path, tumor, 'intermediate_results')
    try:
        t_peps = pd.read_csv(tum_peptides_path, sep = '\t')
    except:
        t_peps = pd.DataFrame(columns = ['isoform', 'isoform_num', 'junc_aa_pos', 'peptides', 'junction', 'frame_selection'])
    t_peps = t_peps.append(full_wt_peptides)
    t_peps = t_peps.reset_index(drop = True)
    t_peps.to_csv('{}/all_peptides.tsv'.format(intermediate_dir), sep = '\t', index = False)

# Create input for pyprsent pipeline
# fasta of isoforms for NetMHCpan input
# peptides spanning the junction to subset NetMHCpan output to peptides of interest
# patient specific junction files to map to NetMHCpan short IDs (this software cannot handle IDs longer than ~15 characters)

output_dir = '{}'.format(args.isoforms_dir)
save_dir = '{}'.format(args.output_directory_fastas)
maf_dir = '{}'.format(args.output_directory_mapping)

patient = '{}'.format(args.patient_sample)
    
# read neopeptides.tsv file
path = os.path.join(output_dir, 'tumors', patient, 'intermediate_results', 'all_peptides.tsv')
df = pd.read_csv(path, sep='\t')
df = df.dropna(subset=['isoform'])
df['junc_aa_pos'] = df['junc_aa_pos'].apply(literal_eval)

# class I
df['class_I_peptides'] = df.apply(lambda x: create_peptides('I', x['isoform'], x['junc_aa_pos']), axis=1)

# class II
df['class_II_peptides'] = df.apply(lambda x: create_peptides('II', x['isoform'], x['junc_aa_pos']), axis=1)

# assign short netmhcpan ID
new_ids = [i for i in range(0, len(df))]
df['NetMHCpan_short_ID'] = new_ids

# save file for mapping later
savepath = os.path.join(maf_dir, '{}.tsv'.format(patient))
df.to_csv(savepath, sep='\t', index=False)

# create fasta/peptides files
fasta_I_savepath = os.path.join(save_dir, 'mhc_i', '{}.fasta'.format(patient))
fasta_II_savepath = os.path.join(save_dir, 'mhc_ii', '{}.fasta'.format(patient))

peptides_I_savepath = os.path.join(save_dir, 'mhc_i', '{}.peptides'.format(patient))    
peptides_II_savepath = os.path.join(save_dir, 'mhc_ii', '{}.peptides'.format(patient))

with open(fasta_I_savepath, 'w') as fasta_I, open(peptides_I_savepath, 'w') as peptide_I,\
open(fasta_II_savepath, 'w') as fasta_II, open(peptides_II_savepath, 'w') as peptide_II:

    for _, row in df.iterrows():
        # mhc I
        if len(row['isoform'])>=8 and len(row['class_I_peptides'])>0:
            fasta_I.write('>{}\n'.format(row['NetMHCpan_short_ID']))
            fasta_I.write('{}\n'.format(row['isoform']))

            peptide_I.write('>{}\n'.format(row['NetMHCpan_short_ID']))
            peptide_I.write('{}\n'.format(','.join(row['class_I_peptides'])))
        # mhc II
        if len(row['isoform'])>=15 and len(row['class_II_peptides'])>0:
            fasta_II.write('>{}\n'.format(row['NetMHCpan_short_ID']))
            fasta_II.write('{}\n'.format(row['isoform']))

            peptide_II.write('>{}\n'.format(row['NetMHCpan_short_ID']))
            peptide_II.write('{}\n'.format(','.join(row['class_II_peptides'])))

# Prepare patient HLA types for NetMHCpan quantification

mhcI_type_df = pd.read_csv('{}'.format(args.allele_table), sep='\t', index_col=0)

todo_i_allele = (','.join(mhcI_type_df.loc[patient].values))

if args.net_two:
    todo_i_allele_list = todo_i_allele.split(',')
    todo_i_allele_list = list(set(todo_i_allele_list))
    todo_i_allele = ','.join(todo_i_allele_list)

# create bash script to run netmhcpan for the given sample
def create_cluster_script(patient, todo_i_allele):

    with open('{}/supplemental/temp_files/{}_netmhcpan.sh'.format(args.pipeline_directory, patient), 'w') as out_file:
        out_file.write("#! /bin/bash\n")
        # Run netmhcpan
        if args.net_two:
            out_file.write('netMHCIIpan -a {} -f {}/mhc_ii/{}.fasta -xls -xlsfile {}/{}.xlsoutput > /dev/null'.format(todo_i_allele, args.output_directory_fastas, patient, args.netmhcpan_output_directory, patient))
        else:
            out_file.write('netMHCpan -a {} -f {}/mhc_i/{}.fasta -xls -xlsfile {}/{}.xlsoutput > /dev/null'.format(todo_i_allele, args.output_directory_fastas, patient, args.netmhcpan_output_directory, patient))
        
create_cluster_script(patient, todo_i_allele)
