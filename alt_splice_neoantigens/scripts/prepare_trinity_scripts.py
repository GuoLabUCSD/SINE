import argparse
import pandas as pd
import os

# get case name and RNA Tumor ID
parser = argparse.ArgumentParser()
path_args = parser.add_argument_group("Input/output options:")

path_args.add_argument('-b', '--bam_filepath', type=str, required=True, 
    help="Path to RNA bam file (aligned to genome)")
path_args.add_argument('-j', '--junction_file', type=str, required=True, 
    help="Path to junctions of interest.") 
path_args.add_argument('-o', '--output_script_path', type=str, required=True, 
    help="Path to output Trinity commands")
path_args.add_argument('-w', '--work_dir', type=str, required=True, 
    help="Directory to place intermediate files")
path_args.add_argument('-T', '--trinity_sif', type=str, required=True, 
    help="Path to Trinity singularity image") 
path_args.add_argument('-P', '--pipeline_dir', type=str, required=True, 
    help="Path to pipeline") 
path_args.add_argument('-r', '--read_boost', type=int, required=False,
    help="Read Duplication Factor")
path_args.add_argument('-c', '--clean_intermediate_files', type=bool, 
    default=True, 
    help="Set to False to retain intermediate generated files")

args = parser.parse_args()

junc_df = pd.read_csv(args.junction_file, sep='\t')


def prepare_cmd(f, junction):   
    f.write('echo " - {}"\n'.format(junction))
    
    f.write(f'mkdir -p {args.work_dir}/{junction}\n')
    f.write(f'cd {args.work_dir}/{junction}\n')

    # first subset to reads spanning the region designated by the junction of interest
    f.write(f'samtools view --fetch-pairs -b {args.bam_filepath} {junction} > {junction}.trinity_in_prefilter.bam\n')
   
    # index the new bam file
    f.write(f'samtools index {junction}.trinity_in_prefilter.bam\n')
 
    # separate out reads with junctions together, just 5', or just 3'
    f.write(f'python {args.pipeline_dir}/scripts/subset_to_junction_reads.py --input_bam_path {junction}.trinity_in_prefilter.bam --output_bam_prefix {junction}.trinity_in --junction {junction}\n')
    
    # sort reads by name for bedtools
    f.write(f'samtools sort -n {junction}.trinity_in.ase.bam > {junction}.trinity_in_sorted.ase.bam\n')

    # convert to fastq
    #f.write(f'bedtools bamtofastq -i {junction}.trinity_in.ase.bam -fq {junction}.trinity_in.ase.fq\n')
    f.write(f'bedtools bamtofastq -i {junction}.trinity_in_sorted.ase.bam -fq {junction}.trinity_in.ase.fq -fq2 {junction}.trinity_in2.ase.fq\n')
    # f.write(f'bedtools bamtofastq -i {junction}.trinity_in.wildtype.bam -fq {junction}.trinity_in.wildtype.fq\n\n')

    if args.read_boost:
        f.write(f'python {args.pipeline_dir}/scripts/read_booster.py --event_dir {args.work_dir}/{junction} --read_boost {args.read_boost}\n')
        f.write(f'if [ -s {junction}.trinity_in_sorted.ase.bam ]; then\n')
        f.write(f'\tsingularity exec -e {args.trinity_sif} Trinity --seqType fq --left {junction}.trinity_in_boosted.ase.fq --right {junction}.trinity_in2_boosted.ase.fq --SS_lib_type RF --max_memory 10G --output trinity_out_{junction}_ase --min_contig_length 50 > /dev/null\n')
        f.write('fi\n')

    else:
        # only run trinity for bams with reads inside
        #f.write(f'if [ -s {junction}.trinity_in.ase.bam ]; then\n')
        #f.write(f'\tsingularity exec -e {args.trinity_sif} Trinity --seqType fq --single {junction}.trinity_in.ase.fq --max_memory 10G --output trinity_out_{junction}_ase --min_contig_length 50 > /dev/null\n')
        f.write(f'if [ -s {junction}.trinity_in_sorted.ase.bam ]; then\n')
        f.write(f'\tsingularity exec -e {args.trinity_sif} Trinity --seqType fq --left {junction}.trinity_in.ase.fq --right {junction}.trinity_in2.ase.fq --SS_lib_type RF --max_memory 10G --output trinity_out_{junction}_ase --min_contig_length 50 > /dev/null\n')
        f.write('fi\n')
    

    ### -- WE DON'T CARE ABOUT WILDTYPE HERE -- ##
    # f.write(f'if [ -s {junction}.trinity_in.wildtype.bam ]; then\n')
    # f.write(f'\tsingularity exec -e {args.trinity_sif} Trinity --seqType fq --single {junction}.trinity_in.wildtype.fq --max_memory 10G --output trinity_out_{junction}_wildtype --min_contig_length 50 > /dev/null\n\n')
    # f.write('fi\n')
    ### -- WE DON'T CARE ABOUT WILDTYPE HERE -- ##
    

    # if bool(args.clean_intermediate_files):
    #     f.write(f'rm {junction}.trinity_in_prefilter.bam\n') 
    #     f.write('rm *fq\n')
    #     f.write(f'rm {junction}.trinity_in.*bam\n\n')
    

    return


# check args are good
assert os.path.isfile(args.bam_filepath)
assert os.path.isfile(args.junction_file)
assert os.path.isfile(args.trinity_sif)

# if work_dir DNE create
if not os.path.isdir(args.work_dir):
    print('\t- Creating {}'.format(args.work_dir))
    os.system('mkdir -p {}'.format(args.work_dir))

# if output path exists already, give warning
if os.path.isfile(args.output_script_path):
    print('\t- WARNING - {} will be overwrriten. Continuing..'.format(args.output_script_path))

# create trinity commands
with open(args.output_script_path, 'w') as f:

    for _, row in junc_df.iterrows():
        prepare_cmd(f, row['junction'])

