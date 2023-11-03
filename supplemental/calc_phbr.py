import argparse
import pandas as pd
import numpy as np
import os
from scipy.stats import hmean as harmonic_mean
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from ast import literal_eval
import re
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

args = parser.parse_args()

# Load junctions we generated neopeptides for
junc_df = pd.read_csv('{}'.format(args.junctions_table), sep='\t')

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

# Prioritization
same_start =  total_df[total_df.duplicated(['chr_intron_start','patient'], keep=False)]
same_start = same_start.sort_values(by=['chr_intron_start','junction_category'], ascending=True)

same_start_total_df = pd.DataFrame()

for patient in same_start['patient'].unique():
    
    pat_df = same_start[same_start['patient']==patient]
    pat_junc_df = pat_df.pivot(index=['chr_intron_start'], columns=['junction_category'], values='junction')
    pat_exp_df = pat_df.pivot_table(index=['chr_intron_start'], columns=['junction_category'], values='num_uniq_reads_at_junc')
    pat_pep_df = pat_df.pivot_table(index=['chr_intron_start'], columns=['junction_category'], values='peptide', aggfunc = np.sum)
    pat_iso_df = pat_df.pivot_table(index=['chr_intron_start'], columns=['junction_category'], values='isoform', aggfunc = np.sum)
    pat_pos_df = pat_df.pivot_table(index=['chr_intron_start'], columns=['junction_category'], values='junc_aa_pos', aggfunc = np.sum)
    pat_df = pat_df.pivot_table(index=['chr_intron_start'], columns=['junction_category'], values='PHBR_scores')

    pat_exp_df = pat_exp_df.rename(columns={'ASE_junc_of_interest': 'ASE_num_uniq', 'other': 'WT_num_uniq'})
    pat_pep_df = pat_pep_df.rename(columns={'ASE_junc_of_interest': 'ASE_peptide', 'other': 'WT_peptide'})
    pat_iso_df = pat_iso_df.rename(columns={'ASE_junc_of_interest': 'ASE_Isoform', 'other': 'WT_Isoform'})
    pat_pos_df = pat_pos_df.rename(columns={'ASE_junc_of_interest': 'ASE_AA_Junc_Position', 'other': 'WT_AA_Junc_Position'})
    pat_df = pat_df.rename(columns={'ASE_junc_of_interest': 'ASE_PHBR', 'other': 'WT_PHBR'})
    
    pat_df['ase_wt_ratio'] = pat_df['ASE_PHBR'] / pat_df['WT_PHBR']
    pat_exp_df['ase_wt_exp_ratio'] = pat_exp_df['ASE_num_uniq'] / pat_exp_df['WT_num_uniq']
    pat_df = pat_df.join(pat_junc_df[['ASE_junc_of_interest']]).join(pat_junc_df[['other']]).join(pat_exp_df).join(pat_pep_df).join(pat_iso_df).join(pat_pos_df)
    pat_df['patient'] = patient
    pat_df = pat_df.rename(columns={'other': 'WT_junc'})
    same_start_total_df = same_start_total_df.append(pat_df)

same_start_total_df['category'] = 'shared_start'

#plt.figure(figsize=(20,20))
#sns.barplot(x='ASE_junc_of_interest', y='ase_wt_ratio', data=same_start_total_df)
#plt.xticks(rotation=45, ha='right')
#plt.savefig('{}/same_start_ratios.png'.format(args.output_directory))
#plt.close()

same_end = total_df[total_df.duplicated(['chr_intron_end','patient'], keep=False)]
same_end_total_df = pd.DataFrame()

for patient in same_end['patient'].unique():
    pat_df = same_end[same_end['patient']==patient]
    pat_junc_df = pat_df.pivot(index=['chr_intron_end'], columns=['junction_category'], values='junction')
    pat_exp_df = pat_df.pivot_table(index=['chr_intron_end'], columns=['junction_category'], values='num_uniq_reads_at_junc')
    pat_pep_df = pat_df.pivot_table(index=['chr_intron_end'], columns=['junction_category'], values='peptide', aggfunc = np.sum)
    pat_iso_df = pat_df.pivot_table(index=['chr_intron_end'], columns=['junction_category'], values='isoform', aggfunc = np.sum)
    pat_pos_df = pat_df.pivot_table(index=['chr_intron_end'], columns=['junction_category'], values='junc_aa_pos', aggfunc = np.sum)
    pat_df = pat_df.pivot_table(index=['chr_intron_end'], columns=['junction_category'], values='PHBR_scores')

    pat_df = pat_df.rename(columns={'ASE_junc_of_interest': 'ASE_PHBR', 'other': 'WT_PHBR'})    
    pat_pep_df = pat_pep_df.rename(columns={'ASE_junc_of_interest': 'ASE_peptide', 'other': 'WT_peptide'})
    pat_exp_df = pat_exp_df.rename(columns={'ASE_junc_of_interest': 'ASE_num_uniq', 'other': 'WT_num_uniq'})
    pat_iso_df = pat_iso_df.rename(columns={'ASE_junc_of_interest': 'ASE_Isoform', 'other': 'WT_Isoform'})
    pat_pos_df = pat_pos_df.rename(columns={'ASE_junc_of_interest': 'ASE_AA_Junc_Position', 'other': 'WT_AA_Junc_Position'})


    pat_df['ase_wt_ratio'] = pat_df['ASE_PHBR'] / pat_df['WT_PHBR']
    pat_exp_df['ase_wt_exp_ratio'] = pat_exp_df['ASE_num_uniq'] / pat_exp_df['WT_num_uniq'] 

    pat_df = pat_df.join(pat_junc_df[['ASE_junc_of_interest']]).join(pat_junc_df[['other']]).join(pat_exp_df).join(pat_pep_df).join(pat_iso_df).join(pat_pos_df)
    pat_df['patient'] = patient
    pat_df = pat_df.rename(columns={'other': 'WT_junc'})
    same_end_total_df = same_end_total_df.append(pat_df)

same_end_total_df['category'] = 'shared_end'

#plt.figure(figsize=(20,20))
#sns.barplot(x='ASE_junc_of_interest', y='ase_wt_ratio', data=same_end_total_df)
#plt.xticks(rotation=45, ha='right')
#plt.savefig('{}/same_end_ratios.png'.format(args.output_directory))
#plt.close()

total_ratio_df = same_start_total_df.append(same_end_total_df)
total_ratio_df = total_ratio_df.sort_values(by=['ase_wt_ratio','ase_wt_exp_ratio', 'ASE_PHBR', 'WT_PHBR'], ascending=[True,False,True,False])

savepath = '{}/fin_prioritization.high_express_junctions.tsv'.format(args.output_directory)
total_ratio_df.to_csv(savepath, sep='\t', index=False)

#x, y = 'ase_wt_ratio','ase_wt_exp_ratio'
#plt.figure(figsize=(10,10))
#sns.scatterplot(x=x, y=y, data=total_ratio_df)

#for _, row in total_ratio_df.iterrows():
#    if row[x]<0.2 and row[y]>2:
#        plt.text(row[x], row[y], '  '+row['ASE_junc_of_interest'])
#plt.yscale('log')
#plt.xscale('log')        
#plt.savefig('{}/scatter_ratios.png'.format(args.output_directory))
#plt.close()

#total_ratio_df = pd.read_csv('{}/fin_prioritization.high_express_junctions.tsv'.format(args.output_directory), sep = '\t')
#plt.figure(figsize=(20,20))
#sns.barplot(data = total_ratio_df, x='ASE_junc_of_interest', y='ase_wt_ratio', hue = 'category')
#plt.xticks(rotation=45, ha='right')
#plt.savefig('{}/start_and_end_ratios.png'.format(args.output_directory))
#plt.close()

# Get event occurance and sort
my_file = pd.read_csv('{}/fin_prioritization.high_express_junctions.tsv'.format(args.output_directory), sep = '\t')
my_file['%patients_with_ASEpeptide']=my_file.groupby(['ASE_peptide'])['ASE_peptide'].transform('count')
my_file['%patients_with_ASEpeptide']=(my_file['%patients_with_ASEpeptide'] / len(patient_list))*100
my_file = my_file.sort_values('ASE_junc_of_interest')
my_file.to_csv('{}/fin_sorted_prioritization.high_express_junctions.tsv'.format(args.output_directory), sep = '\t', index = False)

# Filter Priority File to show junctions where ASE PHBR < WT
my_file = pd.read_csv('{}/fin_prioritization.high_express_junctions.tsv'.format(args.output_directory), sep = '\t')
#Filter just by PHBR
my_file = my_file.loc[my_file['ase_wt_ratio'] < 1]
my_file['%patients_with_ASEpeptide']=my_file.groupby(['ASE_peptide'])['ASE_peptide'].transform('count')
my_file['%patients_with_ASEpeptide']=(my_file['%patients_with_ASEpeptide'] / len(patient_list))*100
my_file = my_file.sort_values('ASE_junc_of_interest')
my_file.to_csv('{}/fin_phbr_filtered_prioritization.high_express_junctions.tsv'.format(args.output_directory), sep='\t', index=False)

#Filter by PHBR and ASE Expression > WT if necessary
my_file = pd.read_csv('{}/fin_prioritization.high_express_junctions.tsv'.format(args.output_directory), sep = '\t')
my_file = my_file.loc[my_file['ase_wt_ratio'] < 1]
my_file = my_file.loc[my_file['ase_wt_exp_ratio'] > 1]
my_file['%patients_with_ASEpeptide']=my_file.groupby(['ASE_peptide'])['ASE_peptide'].transform('count')
my_file['%patients_with_ASEpeptide']=(my_file['%patients_with_ASEpeptide'] / len(patient_list))*100
my_file = my_file.sort_values('ASE_junc_of_interest')
my_file.to_csv('{}/fin_exp_and_phbr_filtered_prioritization.high_express_junctions.tsv'.format(args.output_directory), sep='\t', index=False)

