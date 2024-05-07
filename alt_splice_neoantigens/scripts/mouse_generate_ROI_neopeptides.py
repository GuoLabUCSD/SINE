import argparse
import pandas as pd
import os
import numpy as np
from Bio.SeqIO import parse 
from Bio.Seq import reverse_complement, translate
import sys
from _load_seq_dict import load_CDS_dict
from _load_junction_strand import get_junc_strand
from mouse_slice_retrieve_gtf import *
import pysam
from Bio import pairwise2
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def get_consensus_seqs(c_path, kept_reads, ins_positions, insertion, verbose=False):

    isoform_list = []
    with open(c_path, 'r') as fasta_f:
        for entry in parse(fasta_f, 'fasta'):
            isoform_list.append(str(entry.seq))

    ins_pos_list, consensus_seq_list = [],[]

    for consensus_seq in isoform_list:
        
        # consensus sequence orientation
        if not any([y in consensus_seq for y in kept_reads]):
            consensus_seq = reverse_complement(consensus_seq)
        if verbose:
            print(consensus_seq,'\n')
        if verbose:
            print(insertion)
        
        # get insertion location in consensus sequence
        ins_start, ins_end = np.nan, np.nan
        for i, r in enumerate(kept_reads):
            if r in consensus_seq and ins_positions[i][0] != -1:
                if verbose:
                    print(r, consensus_seq.index(r), ins_positions[i])

                ins_start = consensus_seq.index(r)+ins_positions[i][0]
            if r in consensus_seq and ins_positions[i][1] != -1:

                ins_end = consensus_seq.index(r)+ins_positions[i][1]

        if verbose:
            print(ins_start, ins_end, consensus_seq)
      
        if str(ins_start) == 'nan' and str(ins_end) == 'nan':
            continue
        if str(ins_end) == 'nan':
            start_val = int(insertion.split(':')[1].split('-')[0])-1   
            end_val = int(insertion.split(':')[1].split('-')[1])-1
            ins_length = (end_val - start_val) + 1
            ins_end = (ins_start + ins_length) - 1
        elif str(ins_start) == 'nan':
            start_val = int(insertion.split(':')[1].split('-')[0])-1   
            end_val = int(insertion.split(':')[1].split('-')[1])-1
            ins_length = (end_val - start_val) + 1
            ins_start = (ins_end - ins_length) + 1
 
        ins_pos_list.append((ins_start, ins_end))
        consensus_seq_list.append(consensus_seq)
    return ins_pos_list, consensus_seq_list

def nan_index(l, value):
    '''
    l: list,
    value: item to index
    '''
    for i, item in enumerate(l):
        if value==item:
            return i
    return -1

def get_insertion_read_pos(reads, insertion):
    '''Given pysam AlignedSegments return all reads with insertion and what position it is in the read'''
    start = int(insertion.split(':')[1].split('-')[0])-1      # zero based
    end = int(insertion.split(':')[1].split('-')[1])-1          # zero based
    
    kept_reads = []
    insertion_pos_in_reads = []
    
    for read in reads:
        if start in read.positions or end in read.positions:
            s = reverse_complement(read.get_forward_sequence()) if read.is_reverse else read.get_forward_sequence()

            kept_reads.append(s)
            insertion_pos_in_reads.append((nan_index(read.get_reference_positions(full_length=True), start), nan_index(read.get_reference_positions(full_length=True), end)))
            
    return kept_reads, insertion_pos_in_reads


def create_peptides_from_read(protein_seq, ins_aa_pos, mhc):
    '''Returns entire pep'''
    kmer_list = [8,9,10,11] if mhc == 'I' else [15]
    allowed_indicies = []
    if str(ins_aa_pos[0]) == 'nan' or str(ins_aa_pos[1]) == 'nan':
        return None
    if mhc == 'I':
        for i in range(ins_aa_pos[0], ins_aa_pos[1]+12):
            allowed_indicies.append(i)
    if mhc == 'II':
        for i in range(ins_aa_pos[0], ins_aa_pos[1]+16):
            allowed_indicies.append(i)
    peptides = []
    for kmer in kmer_list:
        for i in range(len(protein_seq) - (kmer-1)):
            start = i
            end = i + kmer
            if end == ins_aa_pos[0]:
                continue
            if end in allowed_indicies and start <= ins_aa_pos[1]:
                peptides.append(protein_seq[start:end])
 
    return peptides
    #else:
    #    return None

def create_peptides_from_consensus(consensus_seq_list, ins_pos_list, insertion, junction='Empty', verbose = False):
    '''
    Depending on insertion location, create peptides
        - start,nan : from start to end of consensus
        - nan, end  : from start of consensus to end
        - start, end: (where start+1=end) +/- buffer spanning insertion
    '''
    total_pep_list = []
    alignment_scores = []
    frame_list = []
    selection_criteria_list = []
    total_aa_list = []
    for i, consensus_seq in enumerate(consensus_seq_list):
        ins_nt_pos = ins_pos_list[i]

        frame, ins_aa, all_alignment_scores, selection_criteria = get_best_frame(consensus_seq, ins_nt_pos, insertion, junction, verbose=False)
        if verbose:
            print(selection_criteria)

        if frame and ins_aa:
            peptides = create_peptides_from_read(frame, ins_aa, mhc='II')
            total_pep_list.append(peptides)
            total_aa_list.append(ins_aa)
        else:
            total_pep_list.append(None)
            total_aa_list.append(None)

        alignment_scores.append(all_alignment_scores)
        frame_list.append(frame)
        selection_criteria_list.append(selection_criteria)
    return total_pep_list, frame_list, total_aa_list, alignment_scores, selection_criteria_list


def get_best_frame(seq, nt_pos, insertion, junction, std_thresh=10, verbose=False):
    '''
    Take best frame if the std is greater than some threshold, if not, take frame with least *
    '''
    translated_seq_frame_list = []
    aa_pos_list = []
    alignment_scores = []
    output_enst_list = []
    
    # check reverse complement too
    rc = reverse_complement(seq)
    rc_nt_pos = ((len(seq)-nt_pos[1])-1,(len(seq)-nt_pos[0])-1)
    
    enst_list = retrieve(slice_gtf(insertion, gtf_path=args.gtf_path))
    if verbose:
        print(enst_list)
    enst_list = [x for x in enst_list if x in CDS_dict]
    enst_list = sorted(enst_list) 
    if verbose:
        print('after filter:',enst_list)

    #Get insertion Strand Info
    
    if args.insertion_file:
        ins_strand = get_junc_strand(args.tumor_rna_bam, args.sample_name, junction)
    else:
        ins_strand = str(0)
    
    for enst in enst_list:
        ref_seq = translate(CDS_dict[enst]) # include stop codons 

        # check strand info, if junction is forward stranded use forward translation frames, reverse stranded used the reverse compliment frames, undefined checks all
        if ins_strand == str(1):
            seqs2do = [(seq, nt_pos)]
        elif ins_strand == str(2):
            seqs2do = [(rc, rc_nt_pos)]
        elif ins_strand == str(0):
            seqs2do = [(rc, rc_nt_pos), (seq, nt_pos)]
        
        for c in seqs2do: 
            temp_seq, temp_nt_pos = c
            
            for i in range(3):

                trans_seq = translate(temp_seq[i:], to_stop=True)

                if len(trans_seq)<6:
                    alignment_scores.append(-1)
                    translated_seq_frame_list.append(trans_seq)
                    aa_pos_list.append((-1,-1))
                    output_enst_list.append(enst)
                    continue


                alignment = pairwise2.align.localms(ref_seq, trans_seq, 2,-5,-1,-0.5, one_alignment_only=True)
                if len(alignment)==0:
                    alignment = [[None,None,-np.inf,None,None]]
                alignment_scores.append(alignment[0][2]/len(trans_seq))
                translated_seq_frame_list.append(trans_seq)

                aa_pos_list.append(((temp_nt_pos[0]-i)//3, (temp_nt_pos[1]-i)//3)) # nan will remain nan
                output_enst_list.append(enst)

    if len(alignment_scores)==0:
        return None, None, None, None
    best_frame = np.argmax(alignment_scores)
    
    id2find = output_enst_list[best_frame]

    if args.insertion_file:
        if strand_dict[id2find] == str(1) and ins_strand == str(1):
            strand_file.write('{}\t{}\tForward_MATCH\n'.format(id2find, junction))
        if strand_dict[id2find] == str(-1) and ins_strand == str(2):
            strand_file.write('{}\t{}\tReverse_MATCH\n'.format(id2find, junction))
        if strand_dict[id2find] == str(1) and ins_strand != str(1):
            strand_file.write('{}\t{}\tMISMATCH\n'.format(id2find, junction))
        if strand_dict[id2find] == str(-1) and ins_strand != str(2):
            strand_file.write('{}\t{}\tMISMATCH\n'.format(id2find, junction))

    gene_file.write('{}\t{}\n'.format(insertion, gene_dict[id2find]))

    return translated_seq_frame_list[best_frame], aa_pos_list[best_frame], alignment_scores, output_enst_list[best_frame]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    path_args = parser.add_argument_group("Input/output options:")

    path_args.add_argument('-i', '--intermediate_results', type=str, required=True, 
        help="Path to intermediate results folder with insertion-specific folders (generated in step 2)")
    path_args.add_argument('-P', '--pipeline_dir', type=str, 
        default="../", 
        help="Path to pipeline") 
    path_args.add_argument('-G', '--gtf_path', type=str, required=True,
        help='Path to GTF file used to grab all ENST IDs overlapping the insertion genomic positions')
    path_args.add_argument('-F', '--cds_fasta_file', type=str, required=True,
        help='Path to the cds FASTA file')
    path_args.add_argument('-s', '--sample_name', type=str, required=True,
        help='Name of the Sample Being Run')
    path_args.add_argument('-b', '--tumor_rna_bam', type=str, required=True,
        help='RNA BAM File Path')
    path_args.add_argument('-j', '--insertion_file', type=str, required=False,
        help='File containing list of insertion regions and {optionally} the corresponding junction region')
   
    args = parser.parse_args()

    # check args are good
    assert os.path.isdir(args.intermediate_results)
    get_CDS_info = load_CDS_dict(args.cds_fasta_file)
    CDS_dict = get_CDS_info[0]
    strand_dict = get_CDS_info[1]
    gene_dict = get_CDS_info[2]

    # final neopeptide output file for all insertions
    total_peptide_df = pd.DataFrame()

    # create a text file where flagging mismatching gene and insertion strand direction
    gene_file = open(args.intermediate_results + '/gene_mappings.txt', 'w')
    gene_file.write('ROI\tSymbol\n')

    if args.insertion_file:
        strand_file = open(args.intermediate_results + '/strand_flags.txt', 'w')
        strand_file.write('CDS_Gene\tSymbol\tROI\tFlag\n')

    for insertion_folder in os.listdir(args.intermediate_results):
        working_dir = os.path.join(args.intermediate_results, insertion_folder)
        if not os.path.isdir(working_dir):
            print(f'\t- skipping {insertion_folder}, does not appear to be insertion output.')
            continue

        insertion = str(insertion_folder)

        if args.insertion_file:
            ins_junc_frame = pd.read_csv('{}'.format(args.insertion_file), sep = '\t')
            junction = ins_junc_frame.loc[ins_junc_frame['ROI'] == insertion, 'junction'].iloc[0]

        ase_bam = os.path.join(working_dir, f'{insertion}.trinity_in_sorted.ase.bam')
        ase_trinity_fasta = os.path.join(working_dir, f'trinity_out_{insertion}_ase', 'Trinity.fasta')
        
        #wt_bam = os.path.join(working_dir, f'{insertion}.trinity_in.wildtype.bam')
        #wt_trinity_fasta = os.path.join(working_dir, f'trinity_out_{insertion}_wildtype','Trinity.fasta')

        if not os.path.isfile(ase_bam):
            print(f'\t-{ase_bam} does not exist. Skipping..')
            continue
        if not os.path.isfile(ase_trinity_fasta):
            print(f'\t-{ase_trinity_fasta} does not exist. Skipping..')
            continue


        
        # main part
        
        if os.stat(ase_bam).st_size == 0: 
            print(f'\t- Skipping {ase_bam} (empty)')
            continue
        # get reads of interest from alignment (for insertion local position)
        bam = pysam.AlignmentFile(ase_bam)
        reads, ins_pos_list = get_insertion_read_pos(bam, insertion)

        # get de novo assembled RNA sequences from Trinity
        if len(reads)>0:

            ins_pos_list, consensus_seq_list = get_consensus_seqs(ase_trinity_fasta, reads, ins_pos_list, insertion)
            if args.insertion_file:
                peptides, frame_list, ins_aa_list, alignment_scores, frame_selection = create_peptides_from_consensus(consensus_seq_list, ins_pos_list, insertion, junction)
            else:
                peptides, frame_list, ins_aa_list, alignment_scores, frame_selection = create_peptides_from_consensus(consensus_seq_list, ins_pos_list, insertion)

        else:
            continue

        # save peptides for each isoform for each insertion
        peptide_df = pd.DataFrame({'isoform': frame_list,
                                   'isoform_num': range(len(consensus_seq_list)),
                                   'ins_aa_pos': ins_aa_list, 
                                   'peptides': peptides, 
                                   'insertion': insertion,
                                   'frame_selection': frame_selection})
        if len(total_peptide_df)==0:
            total_peptide_df = peptide_df.copy()
        else:
            total_peptide_df = total_peptide_df.append(peptide_df)
   
    if args.insertion_file: 
        strand_file.close()
    gene_file.close()
    savepath = os.path.join(args.intermediate_results, 'neopeptides.tsv')
    total_peptide_df.to_csv(savepath, sep='\t', index=False)

