import argparse
import pandas as pd
import numpy as np
import os
from ast import literal_eval
from collections import defaultdict
from load_peptides import load_peptide_dict
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser()
args_path = parser.add_argument_group("Inputs:")

args_path.add_argument('-i', '--netmhcpan_results_directory', type=str, required=True, help='Directory containing xls files with affinity results for all samples from netmhcpan')
args_path.add_argument('-f', '--peptide_directory', type=str, required=True, help='Directory where peptide fasta files are stored from previous step')
args_path.add_argument('-o', '--output_directory', type = str, required=True, help='Location to store phbr score results')
args_path.add_argument('-n', '--net_two', type=str, required=False, help='Use NetMHCIIpan instead')

args = parser.parse_args()

# Parse netmhcpan output
def get_pypresent_output(peptide_filepath, xls_file_list, output_path, affinity_col_base='Rank', output_pep_path=None, verbose=False):
    '''Parses xls output from netMHCpan, saves allele-specific best rank scores
    Args:
        peptide_filepath (str): path to previously saved mutated peptides
        xls_file_list (list): list of paths to NetMHCpan XLS output files
        output_path (str): path to save output file
    '''
    # Load pre-determined mutated peptides
    peptide_dict = load_peptide_dict(peptide_filepath)

    # Prepare allele-specific BR output dataframe
    allele_br_output_df = pd.DataFrame()
    allele_br_pep_output_df = pd.DataFrame()

    # For each XLS file (1 for MHC-I, 3 for MHC-II)
    for xls_output_file in xls_file_list.split(','):
        if verbose:
            print(xls_output_file)
        df = pd.read_csv(xls_output_file, sep='\t', header=1).drop_duplicates()

        # Extract alleles, map to Rank.N column name
        with open(xls_output_file, 'r') as f:
            alleles = f.readline().strip().split()

        allele_rank_dict = defaultdict(str)
        for i, allele in enumerate(alleles): # MHC-II could have duplicates
            if allele in allele_rank_dict.keys():
                continue
            if i == 0:
                allele_rank_dict[allele] = affinity_col_base
            else:
                allele_rank_dict[allele] = '{}.{}'.format(affinity_col_base, i)
        alleles = list(set(alleles))

        # Prepare output dataframe
        mut_rows = []
        kept_muts = []
        br_peps = []
        mutation_IDs = list(set(df.ID.values))
        
        for i, xls_mut_id in enumerate(mutation_IDs):

            # make sure mutation ID in peptide_dict
            if str(xls_mut_id) in peptide_dict.keys():

                # subset to mutation-specific data frame
                temp_mut_df = df[df['ID']==xls_mut_id]

                # subset to allele-specific mutation data frame
                allele_scores = []
                br_pep_list = []
                for allele in alleles:
                    temp_allele_df = temp_mut_df[['Peptide', allele_rank_dict[allele]]]

                    # extract only peptides with residue
                    mutated_peptides = peptide_dict[str(xls_mut_id)]
                    temp_allele_df = temp_allele_df[temp_allele_df['Peptide'].isin(mutated_peptides)]
                    temp_allele_df = temp_allele_df.drop_duplicates()

                    # get scores
                    # try:
                    br_idx = temp_allele_df[allele_rank_dict[allele]].argmin()
                    br_score = temp_allele_df[allele_rank_dict[allele]].values[br_idx]
                    br_pep = temp_allele_df['Peptide'].values[br_idx]
                    allele_scores.append(br_score)
                    br_pep_list.append(br_pep)
                # except:
                #         print('No peptides found for mutation ID: {}, skipping.'.format(xls_mut_id))
                #         pass

                mut_rows.append(allele_scores)
                br_peps.append(br_pep_list)
                kept_muts.append(xls_mut_id)

        output_df = pd.DataFrame(mut_rows, columns=alleles, index=kept_muts)
        output_pep_df = pd.DataFrame(br_peps, columns=alleles, index=kept_muts)

        output_df.index = pd.to_numeric(output_df.index)
        output_pep_df.index = pd.to_numeric(output_pep_df.index)

        if len(allele_br_output_df) == 0:
            allele_br_output_df = output_df.copy()
            allele_br_pep_output_df = output_pep_df.copy()
            #print(allele_br_pep_output_df)
        else:
            allele_br_output_df = pd.concat([allele_br_output_df, output_df], axis=1)
            allele_br_pep_output_df = pd.concat([allele_br_pep_output_df, output_pep_df], axis=1)
    
    allele_br_pep_output_df = allele_br_pep_output_df.add_suffix('_peptide')
    
    allele_br_output_df = allele_br_output_df.join(allele_br_pep_output_df)

    only_scores_df = allele_br_output_df.loc[:,~allele_br_output_df.columns.str.endswith('_peptide')]
    allele_column_names = only_scores_df.columns.to_list()
    only_scores_df['lowest_phbr'] = only_scores_df[allele_column_names].min(axis=1)
    only_scores_df['lowest_allele'] = only_scores_df[allele_column_names].idxmin(axis=1)
    
    allele_br_output_df['lowest_allele'] = only_scores_df['lowest_allele']
    
    allele_br_output_df['lowest_allele'] = allele_br_output_df['lowest_allele'].astype(str) + '_peptide'
    
    lowest_peptide_string = []
    for index, row in allele_br_output_df.iterrows():
        #print(index)
        x = row['lowest_allele']
        lowest_peptide_string.append(row[x])
    allele_br_output_df['peptide'] = lowest_peptide_string
    
    allele_br_output_df = allele_br_output_df[allele_br_output_df.columns.drop(list(allele_br_output_df.filter(regex='_peptide')))]
    allele_br_output_df = allele_br_output_df.drop(['lowest_allele'], axis = 1)    

    #allele_br_output_df.index = allele_br_output_df.index.astype(int)
    #allele_br_output_df.sort_index()
    allele_br_output_df.to_csv(output_path, sep='\t')
    if output_pep_path:
        allele_br_pep_output_df.to_csv(output_pep_path, sep='\t')

xls_output_dir = '{}'.format(args.netmhcpan_results_directory)
peptide_output_dir = '{}'.format(args.peptide_directory)
output_dir = '{}'.format(args.output_directory)

# run for each patient's netmhcpan output file
patient_list = set([x.split('.')[0] for x in os.listdir(xls_output_dir)])
for i, patient in enumerate(patient_list):
    output_path = os.path.join(output_dir, '{}.output'.format(patient))

    peptide_path = os.path.join(peptide_output_dir, '{}.peptides'.format(patient))
    xls_path_list = ','.join([os.path.join(xls_output_dir, x) for x in os.listdir(xls_output_dir) if x.split('.')[0] == patient])

    if args.net_two:
        get_pypresent_output(peptide_path, xls_path_list, output_path, affinity_col_base='Rank')
    else:
        get_pypresent_output(peptide_path, xls_path_list, output_path, affinity_col_base='EL_Rank')

