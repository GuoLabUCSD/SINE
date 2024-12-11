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
    
def keep_insertion_reads_only(bampath, output_prefix, insertion, junction = 'None'):
    '''
    bampath       : path to original bam
    output_prefix : prefix to save the stratified bams to. e.g. <output_dir>/prefix
    insertion      : string following the format: chr1:12344444-12355555     
    '''
    insertion = insertion.split(':')[1]
    start, end = insertion.split('-')
    start = int(start)-1      # zero based
    end = int(end)-1          # zero based
    middle_pos = []
    for i in range(start+1, end):
        middle_pos.append(i)
    used_mates = []
    shared_coord = ''

    #Optionally need to get junction coordinate info and what side of the insertion the junction is connected to
    if junction != 'None':
        junction = junction.split(':')[1]
        j_start, j_end = junction.split('-')
        j_start = int(j_start)-1      # zero based
        j_end = int(j_end)-1          # zero based
        if j_start == end:
            shared_coord = j_start
        if j_end == start:
            shared_coord = j_end
    
    # load bam
    input_bam = pysam.AlignmentFile(bampath, 'rb')
    ase_bam = pysam.AlignmentFile(f'{output_prefix}.ase.bam', 'wb', template=input_bam)

    # check every read
    # write the read mate to the same file and skip the read if the mate has already been written
    for read in input_bam:

        read_header = str(read)

        if read_header in used_mates:
            continue

        if start in read.positions or end in read.positions:
            
            #Optionally need to remove reads that span the region of interest, but not the junction
            if shared_coord == j_end and (start-1) in read.positions:
                continue
            if shared_coord == j_end and (j_start+1) in read.positions:
                continue
            if shared_coord == j_start and (end+1) in read.positions:
                continue
            if shared_coord == j_start and (j_end-1) in read.positions:
                continue
            
            insertion_read_start = nan_index(read.get_reference_positions(full_length=True), start)
            insertion_read_end = nan_index(read.get_reference_positions(full_length=True), end)
            
            # if read mate is unmapped, skip to the next read
            if insertion_read_start != -1 or insertion_read_end != -1:
                try:
                    read_mate = input_bam.mate(read)
                    ase_bam.write(read)
                    ase_bam.write(read_mate)
                    mate_header = str(read_mate)
                    used_mates.append(mate_header)
                except:
                    continue

        elif all(item in middle_pos for item in read.positions):
            try:
                read_mate = input_bam.mate(read)
                ase_bam.write(read)
                ase_bam.write(read_mate)
                mate_header = str(read_mate)
                used_mates.append(mate_header)
            except:
                continue


    # save output bams 
    input_bam.close()
    ase_bam.close()

    return 

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    path_args = parser.add_argument_group("Input options:")

    path_args.add_argument('--input_bam_path', type=str, help='path to input bam', required=True)
    path_args.add_argument('--output_bam_prefix', type=str, help='prefix to save the stratified bams to. e.g. "<output_dir>/some_prefix"', required=True)
    path_args.add_argument("--insertion", type = str, help = "insertion e.g. chr1:111111-111113", required=True)
    path_args.add_argument("--junction", type = str, help = "corresponding junction connected to the insertion", required = False)

    args = parser.parse_args()

    # check inputs
    assert os.path.isfile(args.input_bam_path)

    # warn for overwrite
    if os.path.isfile(args.output_bam_prefix):
        print('WARNING - {} will be overwrriten. Continuing..'.format(args.output_bam_prefix))

    keep_insertion_reads_only(args.input_bam_path, args.output_bam_prefix, args.insertion, args.junction)

