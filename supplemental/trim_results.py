import argparse
import sys
import os
import pandas as pd
from ast import literal_eval
import re

parser = argparse.ArgumentParser()
args_path = parser.add_argument_group("Inputs:")

args_path.add_argument('-o', '--output_directory', type=str, required=True, help='Directory of SINE Results')

args = parser.parse_args()

#Functions
def get_consensus_iso(kmer, my_double_list, my_doubleaas_list):
    fin_list = []
    for val, location in zip(my_double_list, my_doubleaas_list):
        for i, j in zip(val, location):
            j = literal_eval(j)
            if i.find(kmer) >= 0:
                junc_start = j[0]
                junc_end = j[1]
                if junc_start != junc_end:
                    iso_start = i.find(kmer) - 10
                    if iso_start < 0:
                        iso_start = 0
                    kmer_start = i.find(kmer)
                    kmer_end = kmer_start + len(kmer)
                    iso_end = kmer_end + 10
                    i = i[iso_start:junc_start] + '*' + i[junc_start:junc_end+1] + '*' + i[junc_end+1:iso_end]
                    fin_list.append(i)
                else:
                    iso_start = i.find(kmer) - 10
                    if iso_start < 0:
                        iso_start = 0
                    kmer_start = i.find(kmer)
                    kmer_end = kmer_start + len(kmer)
                    iso_end = kmer_end + 10
                    i = i[iso_start:junc_start] + '*' + i[junc_start] + '*' + i[junc_start+1:iso_end]
                    fin_list.append(i)
    iso2keep = max(set(fin_list), key=fin_list.count)
    return iso2keep

def trim_iso(iso):
    firstoccurance = iso.find('*') 
    secondoccurance = iso.find('*', firstoccurance + 1)
    start = iso[0: firstoccurance]
    middle = iso[firstoccurance: secondoccurance+1]
    end = iso[secondoccurance+1::]
    if secondoccurance - firstoccurance == 2:
        if len(start) < 9 or len(end) < 10:
            trim = 'Too short to trim'
        else:
            trim = start[-9::] + middle + end[0:10]
    else:
        if len(start) < 9 or len(end) < 9:
            trim = 'Too short to trim'
        else:
            trim = start[-9::] + middle + end[0:9]
    return trim

#Read files
sine_results_directory = '{}'.format(args.output_directory)

trimmed_df = pd.DataFrame(columns = ['junction', 'Gene', 'Highest Binding Junction Spanning Peptide', 'Consensus Isoform', 'Trimmed 20mer',
                                    'Avg PHBR', 'Avg WT PHBR', '%Samples w/ Junction'])

fin_results = pd.read_csv('{}/fin_results.txt'.format(sine_results_directory), sep = '\t')
raw_results = pd.read_csv('{}/raw_results.txt'.format(sine_results_directory), sep = '\t')

junctionstocheck = fin_results['ASE_junction'].to_list()

count = 0
for i in junctionstocheck:
    fin_row = fin_results.loc[fin_results['ASE_junction'] == i].reset_index(drop = True)
    raw_row = raw_results.loc[raw_results['ASE_junction'] == i].reset_index(drop = True)
    kmer = literal_eval(fin_row['Best_ASE_Peptide_HLA_Pair'][0])[0]
    iso_lists = literal_eval(raw_row['ASE_isoforms'][0])
    index_lists = literal_eval(raw_row['ASE_junction_indices'][0])
    x = get_consensus_iso(kmer, iso_lists, index_lists)
    trim = trim_iso(x)
    trimmed_df.loc[count] = [i, fin_row['Symbol'][0], kmer, x, trim, fin_row['Average_ASE_PHBR_Score(s)'][0], fin_row['Average_WT_PHBR_Score(s)'][0], fin_row['%samples_with_ASE_junction'][0]]
    count = count + 1
trimmed_df.to_csv('{}/trimmed_results.txt'.format(args.output_directory), sep = '\t', index = False)

