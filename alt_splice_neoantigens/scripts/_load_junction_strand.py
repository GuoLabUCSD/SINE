import pandas as pd
import os
import sys

# load CDS/cDNA file

def get_junc_strand(BAM_path, sample_name, junction):
    path = os.path.dirname(BAM_path)
    path = os.path.join(path, sample_name+'SJ.out.tab')
    star_cols = ['chr','intron_start', 'intron_end','strand','intron_motif','in_annot_database','num_uniq_reads_at_junc','num_multi_reads_at_junc','max_spliced_overhang']
    junction_file = pd.read_csv(path, sep = '\t', header = None, low_memory = False)
    junction_file.columns = star_cols
    junction_file['junction'] = junction_file.apply(lambda x: '{}:{}-{}'.format(x['chr'], x['intron_start']-1, x['intron_end']+1), axis=1)
    junction_file = junction_file.drop(columns = ['chr', 'intron_start', 'intron_end', 'intron_motif', 'in_annot_database', 'num_uniq_reads_at_junc', 'num_multi_reads_at_junc', 'max_spliced_overhang'])
    junction_file = junction_file.set_index('junction')
    row = junction_file.loc[junction]
    strand = str(row['strand'])
    return strand

