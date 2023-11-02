import argparse
import pandas as pd
import numpy as np
import os
import re
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')

parser = argparse.ArgumentParser()
args_path = parser.add_argument_group("Inputs:")

args_path.add_argument('-i', '--junctions_list', type=str, required=True, help='Path to a new line separated text file with junction names. The header line should be "junction" and the junctions should be in "#:start-end" format where "#" is the mouse chromosome number, or X/Y.')
args_path.add_argument('-b', '--STAR_directory', type=str, required=True, help='Path to the directory of STAR output (BAM + SJ.out.tab')
args_path.add_argument('-o', '--output_directory', type=str, required=True, help='Location to store output files')
args_path.add_argument('-n', '--output_prefix', type=str, required=True, help='Desired prefix/name for result files')

args = parser.parse_args()

# load data 
junc_df = pd.read_csv('{}'.format(args.junctions_list), sep=',')
junc_df = junc_df[['junction']]
junc_df[['chr', 'range']] = junc_df['junction'].str.split(':', expand = True)
junc_df[['start', 'end']] = junc_df['range'].str.split('-', expand = True)
junc_df = junc_df.astype({"start":"int","end":"int"})
junc_df['start_intron'] = junc_df['start'] + 1
junc_df['end_intron'] = junc_df['end'] - 1
junc_df = junc_df.astype({"start":"str","end":"str", "start_intron":"str", "end_intron":"str"})

# get bams to do
bam_list = os.listdir('{}'.format(args.STAR_directory))
bam_list = [x for x in bam_list if x.endswith('.sortedByCoord.out.bam')]

# get patient list
patient_list = [x.split('Aligned')[0] for x in bam_list]

### Look for other junctions that share either a start or end with our junction of interest
cols_of_interest = ['junction','chr','start','end','start_intron','end_intron']
junc_df = junc_df[cols_of_interest]

# create chr-start and chr-end columns
junc_df['chr_intron_start'] = junc_df.apply(lambda x: '{}_{}'.format(x['chr'], x['start_intron']), axis=1)
junc_df['chr_intron_end'] = junc_df.apply(lambda x: '{}_{}'.format(x['chr'], x['end_intron']), axis=1)

# load SJ.out.tab files from Mice and identify any WT junctions w/ a shared start|end 
# SJ.out.tab files are missing "chr" in front of the chromosome number in the first column, edit that here
# SJ.out.tab files also contain un-needed junctions to be removed
star_cols = ['chr','intron_start', 'intron_end','strand','intron_motif','in_annot_database','num_uniq_reads_at_junc','num_multi_reads_at_junc','max_spliced_overhang']

output_df = pd.DataFrame()
for patient in patient_list:
    path = '{}/{}SJ.out.tab'.format(args.STAR_directory, patient)
    star_df = pd.read_csv(path, sep='\t', header=None, low_memory = False)
    star_df.columns = star_cols
    star_df = star_df[star_df["chr"].str.contains("JH") == False]
    star_df = star_df[star_df["chr"].str.contains("GL") == False]
    star_df = star_df[star_df["chr"].str.contains("MU") == False]
    
    # identify shared start/end
    star_df['chr_intron_start'] = star_df.apply(lambda x: '{}_{}'.format(x['chr'], x['intron_start']), axis=1)
    star_df['chr_intron_end'] = star_df.apply(lambda x: '{}_{}'.format(x['chr'], x['intron_end']), axis=1)

    # subset rows to keep 
    star_df = star_df[(star_df['chr_intron_start'].isin(junc_df['chr_intron_start']))|(star_df['chr_intron_end'].isin(junc_df['chr_intron_end']))]
    if len(star_df)>0:
        star_df['junction'] = star_df.apply(lambda x: '{}:{}-{}'.format(x['chr'], x['intron_start']-1, x['intron_end']+1), axis=1)
        found_junctions = star_df[star_df['junction'].isin(junc_df['junction'].values)]['junction'].values
        star_df['junction_category'] = star_df['junction'].apply(lambda x: 'ASE_junc_of_interest' if x in found_junctions  else 'other')
        star_df['patient'] = patient
        output_df = output_df.append(star_df)

# keep highest WT junction
output_df = output_df.sort_values(by=['num_uniq_reads_at_junc'], ascending=False)
output_df = output_df.drop_duplicates(subset=['chr_intron_start','patient','junction_category'], keep='first')
output_df = output_df.drop_duplicates(subset=['chr_intron_end','patient','junction_category'], keep='first')

per_value = np.percentile(output_df[output_df['junction_category']=='other']['num_uniq_reads_at_junc'].values, 75)

# keep highest expressed junctions
output_df = output_df[output_df['num_uniq_reads_at_junc']>per_value]

# save df
savepath = '{}/{}_junctionsOfInterest_plusWT.tsv'.format(args.output_directory, args.output_prefix)
output_df.to_csv(savepath, sep='\t', index=False)

# create junction list per sample
my_junctions = pd.read_csv("{}/{}_junctionsOfInterest_plusWT.tsv".format(args.output_directory, args.output_prefix), sep = '\t')
my_junctions = my_junctions.drop(['chr', 'intron_start', 'intron_end', 'strand', 'intron_motif', 'in_annot_database', 'num_uniq_reads_at_junc', 'num_multi_reads_at_junc', 'max_spliced_overhang', 'chr_intron_start', 'chr_intron_end', 'junction_category'], axis = 1)
patients = my_junctions.groupby('patient')
for name, group in patients:
    group = group.drop(['patient'], axis = 1)
    group.to_csv('{}/{}_junctions.tsv'.format(args.output_directory, name), sep = '\t', index = False)
