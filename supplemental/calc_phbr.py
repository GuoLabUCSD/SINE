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
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
from statannot import add_stat_annotation
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
from plot_medians_in_boxplot import plot_medians_in_boxplot
from stat_annot_utils import create_pairs, create_hue_pairs

parser = argparse.ArgumentParser()
args_path = parser.add_argument_group("Inputs:")

args_path.add_argument('-j', '--junctions_table', type=str, required=True, help='Full junctions table generated from get_wt.sh')
args_path.add_argument('-a', '--affinity_scores', type=str, required=True, help='Directory with allele affinity scores. Ie) Output files from xls_parser.sh')
args_path.add_argument('-m', '--mapping_directory', type=str, required=True, help='Directory containing patient files with the mapping ID number to netmhcpan')
args_path.add_argument('-o', '--output_directory', type=str, required=True, help='Location to store phbr score results and plots')
args_path.add_argument('-t', '--paired_table', type=str, required=True, help='Location of the ASE junctions table with their corresponding paired WT junctions')

args = parser.parse_args()

# Load junctions we generated neopeptides for
junc_df = pd.read_csv('{}'.format(args.junctions_table), sep='\t')
samples_only_df = junc_df.drop_duplicates(subset = 'patient')
tums = samples_only_df[samples_only_df['junction_category'] == 'ASE_junc_of_interest']
norms = samples_only_df[samples_only_df['junction_category'] == 'other']

# Gather NETMHCpan parsed output affinity files for each patient
patient_list = junc_df['patient'].unique()

total_df = pd.DataFrame()
for patient in patient_list:
    aff_path = '{}/{}.output'.format(args.affinity_scores, patient)
    maf_path = '{}/{}.tsv'.format(args.mapping_directory, patient)

    # load file to map NetMHCpan short IDs to the actual junction names and isoform
    pat_df = pd.read_csv(maf_path, sep='\t', usecols=['junction','isoform','NetMHCpan_short_ID', 'junc_aa_pos'])

    # load patient affinity file
    pat_aff = pd.read_csv(aff_path, sep='\t', index_col=0)

    # join
    pat_aff = pat_aff.join(pat_df).set_index('NetMHCpan_short_ID')

    #Drop Peptide Columns to assign back later
    peptide_df = pat_aff[['junction','isoform', 'peptide', 'junc_aa_pos']]
    pat_aff = pat_aff.drop(['peptide', 'isoform', 'junc_aa_pos', 'junction'], axis=1)

    # take harmonic mean as in Marty et al. 2017
    pat_aff = pat_aff.apply(harmonic_mean, axis=1)
    pat_aff = pd.DataFrame(pat_aff, columns=['PHBR_scores']).reset_index()

    #Join Peptides back to the df
    peptide_df = peptide_df.reset_index()
    #print(peptide_df)
    pat_aff = pd.merge(pat_aff, peptide_df, on = 'NetMHCpan_short_ID', how = 'left')

    clean_pat_aff = pd.DataFrame(columns = ['NetMHCpan_short_ID', 'PHBR_scores', 'peptide', 'isoform', 'junc_aa_pos', 'junction'])
    new_index = 0
    for index, row in pat_aff.iterrows():
        clean_pos = literal_eval(row['junc_aa_pos'])
        result = re.finditer(r'(?=({}))'.format(row['peptide']), row['isoform'])
        for match in result:
            if clean_pos[0] in range(match.start(1), match.end(1)) and clean_pos[1] in range(match.start(1), match.end(1)):
                clean_pat_aff.loc[new_index] = row
                new_index = new_index + 1
                continue

    pat_junc_df = junc_df[junc_df['patient']==patient]

    # add manual annotations
    clean_pat_aff = pd.merge(clean_pat_aff, pat_junc_df[['junction_category','junction','chr_intron_start','chr_intron_end','num_uniq_reads_at_junc','patient']],
                       on='junction', how='left')
    clean_pat_aff = clean_pat_aff.drop(['NetMHCpan_short_ID'], axis = 1)

    clean_pat_aff['PHBR_scores'] = pd.to_numeric(clean_pat_aff['PHBR_scores'])


    clean_pat_aff = clean_pat_aff[clean_pat_aff['PHBR_scores'] == clean_pat_aff.groupby('junction')['PHBR_scores'].transform(min)]


    clean_pat_aff = clean_pat_aff.groupby(['PHBR_scores', 'peptide', 'junction', 'junction_category', 'chr_intron_start', 'chr_intron_end', 
                                           'num_uniq_reads_at_junc', 'patient'])[['isoform', 'junc_aa_pos']].agg(list).reset_index()

    total_df = total_df.append(clean_pat_aff)

def plot_boxplots(df, x, y=None, y_list=None, order=None, ylim_dict=None, figsize=(7,4), hue=None, palette='Set1',
        savepath=None, suptitle=None, comparisons_correction=None, plot_significant_only=True):
        if y and y_list:
                print('ERROR. Cannot have y and y_list at the same time')
                return


        plt.figure(figsize=figsize)
        if order is None:
                order = df[x].unique()


        if y_list:

                i=1

                for y in y_list:
                        plt.subplot(1,len(y_list),i)
                        ax = sns.boxplot(x=x, y=y, data=df, order=order, palette=palette)

                        if ylim_dict:
                                if y in ylim_dict:
                                        plt.ylim(ylim_dict[y])

                        add_stat_annotation(ax=ax, x=x, y=y, data=df, test='Mann-Whitney', text_format='simple', order=order, plot_significant_only=plot_significant_only,
                                                                box_pairs=create_pairs(order), verbose=0, loc='outside', comparisons_correction=comparisons_correction)
                        plot_medians_in_boxplot(ax=ax, df=df, x=x, y=y)
                        counts = df[x].value_counts()
                        ax.set_xticklabels(['{}\n({})'.format(x.get_text(),counts.loc[x.get_text()]) for x in ax.get_xticklabels()])
                        plt.xlabel('')
                        i+=1
                plt.tight_layout()

                if suptitle:
                        plt.suptitle(suptitle, y=1.05)

        else:
                if hue:
                        hue_order = df[hue].unique()

                        ax = sns.boxplot(x=x, y=y, data=df, order=order, hue=hue, hue_order=hue_order, palette=palette)

                        if ylim_dict:
                                if y in ylim_dict:
                                        plt.ylim(ylim_dict[y])

                        add_stat_annotation(ax=ax, x=x, y=y, data=df, test='Mann-Whitney', text_format='simple', order=order, hue=hue, hue_order=hue_order,
                                                                box_pairs=create_hue_pairs(order, hue_order), verbose=0, loc='outside', comparisons_correction=comparisons_correction)
                        if len(hue_order)>2:
                                print('Can only handle plotting medians in hue pairs of 2.')
                        else:
                                plot_medians_in_boxplot(ax=ax, df=df, x=x, y=y, hue_left_right=hue_order, hue=hue)
                                counts = df[[x,hue]].value_counts()
                                ax.set_xticklabels([f'{x.get_text()}\n({counts.loc[x.get_text()][hue_order[0]]},{counts.loc[x.get_text()][hue_order[1]]})' for x in ax.get_xticklabels()])
                        plt.xlabel('')
                else:
                        ax = sns.boxplot(x=x, y=y, data=df, order=order, palette=palette)

                        if ylim_dict:
                                if y in ylim_dict:
                                        plt.ylim(ylim_dict[y])

                        add_stat_annotation(ax=ax, x=x, y=y, data=df, test='Mann-Whitney', text_format='simple', order=order,
                                                                box_pairs=create_pairs(order), verbose=0, loc='outside', comparisons_correction=comparisons_correction)
                        plot_medians_in_boxplot(ax=ax, df=df, x=x, y=y)
                        counts = df[x].value_counts()
                        ax.set_xticklabels(['{}\n({})'.format(x.get_text(),counts.loc[x.get_text()]) for x in ax.get_xticklabels()])
                        plt.xlabel('')
        if savepath:
                plt.savefig(savepath, bbox_inches='tight')
                plt.close()
        return

plot_boxplots(x='junction_category', y='PHBR_scores', df=total_df.dropna(), figsize=(3,3), ylim_dict={'PHBR_scores':(-1,21)}, 
              order=['other','ASE_junc_of_interest'], savepath='{}/ase_vs_other.png'.format(args.output_directory))

#Split the combined dataframe of ASE and WT into separate frames
ase_frame = total_df[total_df['junction_category'] == 'ASE_junc_of_interest']
ase_frame = ase_frame.drop(columns = ['junction_category'], axis = 1)
wt_frame = total_df[total_df['junction_category'] == 'other']
wt_frame = wt_frame.drop(columns = ['junction_category'], axis = 1)

all_ases = ase_frame.groupby('junction')
all_wts = wt_frame.groupby('junction')

#Group both dataframes so there is just one junction per row and each cell contains a list of all values found in each patient
final_ase_df = pd.DataFrame(columns = ['ASE_PHBR_Score', 'ASE_peptide', 'ASE_junction', 'ASE_#Reads', 'Tumor_samples', '%samples_with_ASE_junction', 'ASE_isoforms', 'ASE_junction_index'])
final_wt_df = pd.DataFrame(columns = ['WT_PHBR_Score', 'WT_peptide', 'WT_junction', 'WT_#Reads', 'Normal_samples', '%samples_with_WT_junction', 'WT_isoforms', 'WT_junction_index'])

row_index = 0
for name, group in all_ases:
    scores = group['PHBR_scores'].tolist()
    peptides = group['peptide'].tolist()
    junctions = name
    reads = group['num_uniq_reads_at_junc'].to_list()
    samples = group['patient'].to_list()
    all_isos = group['isoform'].to_list()
    perc = (len(samples) / len(tums)) * 100
    all_positions = group['junc_aa_pos'].to_list()
    final_ase_df.loc[row_index] = [scores, peptides, junctions, reads, samples, perc, all_isos, all_positions]
    row_index = row_index + 1

row_index = 0
for name, group in all_wts:
    scores = group['PHBR_scores'].tolist()
    peptides = group['peptide'].tolist()
    junctions = name
    reads = group['num_uniq_reads_at_junc'].to_list()
    samples = group['patient'].to_list()
    all_isos = group['isoform'].to_list()
    perc = (len(samples) / len(norms)) * 100
    all_positions = group['junc_aa_pos'].to_list()
    final_wt_df.loc[row_index] = [scores, peptides, junctions, reads, samples, perc, all_isos, all_positions]
    row_index = row_index + 1

final_ase_df = final_ase_df.set_index('ASE_junction', drop = False)
final_wt_df = final_wt_df.set_index('WT_junction', drop = False)

#Load Pairing Table
pairing_table = pd.read_csv('{}'.format(args.paired_table), sep = '\t')
pairing_table['wt'] = pairing_table['wt'].fillna(0)

#For every row in the pairing table, merge the corresponding ase row with the corresponding wt row
my_final_df = pd.DataFrame(columns = ['ASE_junction', 'ASE_peptide', 'ASE_PHBR_Score', 'ASE_#Reads', 'Tumor_samples', '%samples_with_ASE_junction', 'ASE_isoforms', 'ASE_junction_indices',
                                      'WT_junction', 'WT_peptide', 'WT_PHBR_Score', 'WT_#Reads', 'Normal_samples', '%samples_with_WT_junction', 'WT_isoforms', 'WT_junction_indices'])

index_val = 0
for index, row in pairing_table.iterrows():
    if row['ase'] in final_ase_df.index:
        if row['wt'] == 0 or row['wt'] not in final_wt_df.index:
            x = final_ase_df.loc[row['ase']]
            val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12, val13, val14, val15, val16 = x['ASE_junction'], x['ASE_peptide'], x['ASE_PHBR_Score'], x['ASE_#Reads'], x['Tumor_samples'], x['%samples_with_ASE_junction'], x['ASE_isoforms'], x['ASE_junction_index'], 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan'
            my_final_df.loc[index_val] = val1,val2,val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,val15,val16
            index_val = index_val + 1
        else:
            x = final_ase_df.loc[row['ase']]
            y = final_wt_df.loc[row['wt']]
            val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12, val13, val14, val15, val16 = x['ASE_junction'], x['ASE_peptide'], x['ASE_PHBR_Score'], x['ASE_#Reads'], x['Tumor_samples'], x['%samples_with_ASE_junction'], x['ASE_isoforms'], x['ASE_junction_index'], y['WT_junction'], y['WT_peptide'], y['WT_PHBR_Score'], y['WT_#Reads'], y['Normal_samples'], y['%samples_with_WT_junction'], y['WT_isoforms'], y['WT_junction_index']
            my_final_df.loc[index_val] = val1,val2,val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,val15,val16
            index_val = index_val + 1
    else:
        continue

#Same as above, but make a summarized output file
summary_ase_df = pd.DataFrame(columns = ['ASE_PHBR_Score', 'ASE_peptide', 'ASE_junction', 'ASE_Median#Reads', 'Tumor_samples', '%samples_with_ASE_junction', 'ASE_isoforms'])
summary_wt_df = pd.DataFrame(columns = ['WT_PHBR_Score', 'WT_peptide', 'WT_junction', 'WT_Median#Reads', 'Normal_samples', '%samples_with_WT_junction', 'WT_isoforms'])

row_index = 0
for name, group in all_ases:
    scores = list(set(group['PHBR_scores'].tolist()))
    peptides = list(set(group['peptide'].tolist()))
    junctions = name
    reads = group['num_uniq_reads_at_junc'].to_list()
    median_r = statistics.median(reads)
    samples = group['patient'].to_list()
    isos = group['isoform']
    clean_isos = []
    for i in isos:
        for j in i:
            clean_isos.append(j)
    isos = list(set(clean_isos))
    perc = (len(samples) / len(tums)) * 100
    summary_ase_df.loc[row_index] = [scores, peptides, junctions, median_r, samples, perc, isos]
    row_index = row_index + 1

for name, group in all_wts:
    scores = list(set(group['PHBR_scores'].tolist()))
    peptides = list(set(group['peptide'].tolist()))
    junctions = name
    reads = group['num_uniq_reads_at_junc'].to_list()
    median_r = statistics.median(reads)            
    samples = group['patient'].to_list()
    isos = group['isoform']
    clean_isos = []
    for i in isos:
        for j in i:
            clean_isos.append(j)
    isos = list(set(clean_isos))
    perc = (len(samples) / len(tums)) * 100
    all_positions = group['junc_aa_pos'].to_list()
    summary_wt_df.loc[row_index] = [scores, peptides, junctions, median_r, samples, perc, isos]
    row_index = row_index + 1
summary_ase_df = summary_ase_df.set_index('ASE_junction', drop = False)
summary_wt_df = summary_wt_df.set_index('WT_junction', drop = False)

#Same as above, but for the summarized table
my_summarized_df = pd.DataFrame(columns = ['ASE_junction', 'ASE_peptide', 'ASE_PHBR_Score', 'ASE_Median#Reads', 'Tumor_samples', '%samples_with_ASE_junction', 'ASE_isoforms',
                                      'WT_junction', 'WT_peptide', 'WT_PHBR_Score', 'WT_Median#Reads', 'Normal_samples', '%samples_with_WT_junction', 'WT_isoforms'])

index_val = 0
for index, row in pairing_table.iterrows():
    if row['ase'] in summary_ase_df.index:
        if row['wt'] == 0 or row['wt'] not in summary_wt_df.index:
            x = summary_ase_df.loc[row['ase']]
            val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12, val13, val14 = x['ASE_junction'], x['ASE_peptide'], x['ASE_PHBR_Score'], x['ASE_Median#Reads'], x['Tumor_samples'], x['%samples_with_ASE_junction'], x['ASE_isoforms'], 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan'
            my_summarized_df.loc[index_val] = val1,val2,val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14
            index_val = index_val + 1
        else:
            x = summary_ase_df.loc[row['ase']]
            y = summary_wt_df.loc[row['wt']]
            val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12, val13, val14 = x['ASE_junction'], x['ASE_peptide'], x['ASE_PHBR_Score'], x['ASE_Median#Reads'], x['Tumor_samples'], x['%samples_with_ASE_junction'], x['ASE_isoforms'], y['WT_junction'], y['WT_peptide'], y['WT_PHBR_Score'], y['WT_Median#Reads'], y['Normal_samples'], y['%samples_with_WT_junction'], y['WT_isoforms']
            my_summarized_df.loc[index_val] = val1,val2,val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14
            index_val = index_val + 1
    else:
        continue

#Save Tables
my_final_df.to_csv('{}/fin_results.tsv'.format(args.output_directory), sep='\t', index=False)
my_summarized_df.to_csv('{}/fin_summarized_results.tsv'.format(args.output_directory), sep = '\t', index=False)

