import pysam
import argparse
import os

def nan_index(l, value):
    ''' Helper function to index a list <l> with <value> without returning a NaN
    l: list,
    value: item to index
    '''
    for i, item in enumerate(l):
        if value==item:
            return i
    return -1
    
def keep_junction_reads_only(bampath, output_prefix, junction):
    '''
    bampath       : path to original bam
    output_prefix : prefix to save the stratified bams to. e.g. <output_dir>/prefix
    junction      : string following the format: chr1:12344444-12355555     
    '''
    junction = junction.split(':')[1]
    start, end = junction.split('-')
    start = int(start)-1      # zero based
    end = int(end)-1          # zero based
    
    # load bam
    input_bam = pysam.AlignmentFile(bampath, 'rb')
    ase_bam = pysam.AlignmentFile(f'{output_prefix}.ase.bam', 'wb', template=input_bam)

    # check every read
    for read in input_bam:

        if start in read.positions and end in read.positions:
            
            junc_read_start = nan_index(read.positions, start)
            junc_read_end = nan_index(read.positions, end)
            
            # check if read spans AS junction
            if junc_read_start != -1 and junc_read_end != -1 and junc_read_start+1 == junc_read_end:
                ase_bam.write(read)


    # save output bams 
    input_bam.close()
    ase_bam.close()

    return 

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    path_args = parser.add_argument_group("Input options:")

    path_args.add_argument('--input_bam_path', type=str, help='path to input bam', required=True)
    path_args.add_argument('--output_bam_prefix', type=str, help='prefix to save the stratified bams to. e.g. "<output_dir>/some_prefix"', required=True)
    path_args.add_argument("--junction", type = str, help = "junction e.g. chr1:111111-111113", required=True)

    args = parser.parse_args()

    # check inputs
    assert os.path.isfile(args.input_bam_path)

    # warn for overwrite
    if os.path.isfile(args.output_bam_prefix):
        print('WARNING - {} will be overwrriten. Continuing..'.format(args.output_bam_prefix))

    keep_junction_reads_only(args.input_bam_path, args.output_bam_prefix, args.junction)

