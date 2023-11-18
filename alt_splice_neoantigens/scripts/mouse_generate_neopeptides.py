import argparse
import pandas as pd
import os
import numpy as np
from Bio.SeqIO import parse 
from Bio.Seq import reverse_complement, translate
import sys
from _load_seq_dict import load_CDS_dict
from mouse_slice_retrieve_gtf import *
import pysam
from Bio import pairwise2 

def get_consensus_seqs(c_path, kept_reads, junc_positions, junction, verbose=False):

    isoform_list = []
    with open(c_path, 'r') as fasta_f:
        for entry in parse(fasta_f, 'fasta'):
            isoform_list.append(str(entry.seq))

    junc_pos_list, consensus_seq_list = [],[]

    for consensus_seq in isoform_list:
        
        # consensus sequence orientation
        if not any([y in consensus_seq for y in kept_reads]):
            consensus_seq = reverse_complement(consensus_seq)
        if not any([y in consensus_seq for y in kept_reads]): # if no reads in consensus seq skip
            continue
        if verbose:
            print(consensus_seq,'\n')


        if verbose:
            print(junction)
        
        # get junction location in consensus sequence
        junc_start, junc_end = np.nan, np.nan
        for i, r in enumerate(kept_reads):
            if r in consensus_seq and junc_positions[i][0] != -1:
                if verbose:
                    print(r, consensus_seq.index(r), junc_positions[i])

                junc_start = consensus_seq.index(r)+junc_positions[i][0]
            if r in consensus_seq and junc_positions[i][1] != -1:

                junc_end = consensus_seq.index(r)+junc_positions[i][1]

        if verbose:
            print(junc_start, junc_end, consensus_seq)
        
        junc_pos_list.append((junc_start, junc_end))
        consensus_seq_list.append(consensus_seq)
        
    return junc_pos_list, consensus_seq_list



def nan_index(l, value):
    '''
    l: list,
    value: item to index
    '''
    for i, item in enumerate(l):
        if value==item:
            return i
    return -1
def get_junction_read_pos(reads, junction):
    '''Given pysam AlignedSegments return all reads with junction and what position it is in the read'''
    start = int(junction.split(':')[1].split('-')[0])-1      # zero based
    end = int(junction.split(':')[1].split('-')[1])-1          # zero based
    
    kept_reads = []
    junction_pos_in_reads = []
    
    for read in reads:
        if start in read.positions and end in read.positions:
            s = reverse_complement(read.get_forward_sequence()) if read.is_reverse else read.get_forward_sequence()

            kept_reads.append(s)
            junction_pos_in_reads.append((nan_index(read.positions, start), nan_index(read.positions, end)))
            
    return kept_reads, junction_pos_in_reads


def create_peptides_from_read(protein_seq, junc_aa_pos, mhc):
    '''Returns entire pep'''
    kmer_list = [8,9,10,11] if mhc == 'I' else [15]
    peptides = []
    
    if pd.isna(junc_aa_pos[0]) or pd.isna(junc_aa_pos[1]):
        return None
    
    # make sure junction of interest
    if junc_aa_pos[0]+1 == junc_aa_pos[1] or junc_aa_pos[0]==junc_aa_pos[1]:
        
        for kmer in kmer_list:
            for i in range(len(protein_seq) - (kmer-1)):
                start = i
                end = i + kmer 
            
                if start <= junc_aa_pos[0] and junc_aa_pos[1] < end:
                    peptides.append(protein_seq[start:end])


        return peptides
    else:
        return None

def create_peptides_from_consensus(consensus_seq_list, junc_pos_list, junction, verbose=False):
    '''
    Depending on junction location, create peptides
        - start,nan : from start to end of consensus
        - nan, end  : from start of consensus to end
        - start, end: (where start+1=end) +/- buffer spanning junction
    '''
    total_pep_list = []
    alignment_scores = []
    frame_list = []
    selection_criteria_list = []
    for i, consensus_seq in enumerate(consensus_seq_list):
        junc_nt_pos = junc_pos_list[i]

        frame, junc_aa, all_alignment_scores, selection_criteria = get_best_frame(consensus_seq, junc_nt_pos, junction)
        if verbose:
            print(selection_criteria)

        if frame and junc_aa:
            peptides = create_peptides_from_read(frame, junc_aa, mhc='II')
            total_pep_list.append(peptides)
        else:
            total_pep_list.append(None)
        alignment_scores.append(all_alignment_scores)
        frame_list.append(frame)
        selection_criteria_list.append(selection_criteria)
    return total_pep_list, frame_list, alignment_scores, selection_criteria_list


def create_peptides_from_consensus(consensus_seq_list, junc_pos_list, junction, verbose=False):
    '''
    Depending on junction location, create peptides
        - start,nan : from start to end of consensus
        - nan, end  : from start of consensus to end
        - start, end: (where start+1=end) +/- buffer spanning junction
    '''
    total_pep_list = []
    alignment_scores = []
    frame_list = []
    selection_criteria_list = []
    total_aa_list = []
    for i, consensus_seq in enumerate(consensus_seq_list):
        junc_nt_pos = junc_pos_list[i]

        frame, junc_aa, all_alignment_scores, selection_criteria = get_best_frame(consensus_seq, junc_nt_pos, junction, verbose=False)
        if verbose:
            print(selection_criteria)

        if frame and junc_aa:
            peptides = create_peptides_from_read(frame, junc_aa, mhc='II')
            total_pep_list.append(peptides)
            total_aa_list.append(junc_aa)
        else:
            total_pep_list.append(None)
            total_aa_list.append(None)

        alignment_scores.append(all_alignment_scores)
        frame_list.append(frame)
        selection_criteria_list.append(selection_criteria)
    return total_pep_list, frame_list, total_aa_list, alignment_scores, selection_criteria_list


def get_best_frame(seq, nt_pos, junction, std_thresh=10, verbose=False):
    '''
    Take best frame if the std is greater than some threshold, if not, take frame with least *
    '''
    translated_seq_frame_list = []
    aa_pos_list = []
    alignment_scores = []
    output_enst_list = []
    
    # check reverse complement too
    rc = reverse_complement(seq)
    rc_nt_pos = (len(seq)-nt_pos[1],len(seq)-nt_pos[0])
    
    enst_list = retrieve(slice_gtf(junction, gtf_path=args.gtf_path))
    #print(enst_list)
    if verbose:
        print(enst_list)
    enst_list = [x for x in enst_list if x in CDS_dict]
    #print(enst_list)
    if verbose:
        print('after filter:',enst_list)    
    for enst in enst_list:
        ref_seq = translate(CDS_dict[enst]) # include stop codons 
        
        # for all 6 frames check best transcript alignment
        for c in [(rc, rc_nt_pos), (seq, nt_pos)]: 
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

    return translated_seq_frame_list[best_frame], aa_pos_list[best_frame], alignment_scores, output_enst_list[best_frame]






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    path_args = parser.add_argument_group("Input/output options:")

    path_args.add_argument('-i', '--intermediate_results', type=str, required=True, 
        help="Path to intermediate results folder with junction-specific folders (generated in step 2)")
    path_args.add_argument('-P', '--pipeline_dir', type=str, 
        default="../", 
        help="Path to pipeline") 
    path_args.add_argument('-G', '--gtf_path', type=str, required=True,
        help='Path to GTF file used to grab all ENST IDs overlapping the junction genomic positions')
    path_args.add_argument('-F', '--cds_fasta_file', type=str, required=True,
        help='Path to the cds FASTA file')
   
    args = parser.parse_args()

    # check args are good
    assert os.path.isdir(args.intermediate_results)
    CDS_dict = load_CDS_dict(args.cds_fasta_file)    

    # final neopeptide output file for all junctions
    total_peptide_df = pd.DataFrame()


    for junction_folder in os.listdir(args.intermediate_results):
        working_dir = os.path.join(args.intermediate_results, junction_folder)
        if not os.path.isdir(working_dir):
            print(f'\t- skipping {junction_folder}, does not appear to be junction output.')
            continue

        junction = str(junction_folder)
        ase_bam = os.path.join(working_dir, f'{junction}.trinity_in.ase.bam')
        ase_trinity_fasta = os.path.join(working_dir, f'trinity_out_{junction}_ase', 'Trinity.fasta')
        
        wt_bam = os.path.join(working_dir, f'{junction}.trinity_in.wildtype.bam')
        wt_trinity_fasta = os.path.join(working_dir, f'trinity_out_{junction}_wildtype','Trinity.fasta')

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
        # get reads of interest from alignment (for junction local position)
        bam = pysam.AlignmentFile(ase_bam)
        reads, junc_pos_list = get_junction_read_pos(bam, junction)

        # get de novo assembled RNA sequences from Trinity
        if len(reads)>0:

            junc_pos_list, consensus_seq_list = get_consensus_seqs(ase_trinity_fasta, reads, junc_pos_list, junction)
            peptides, frame_list, junc_aa_list, alignment_scores, frame_selection = create_peptides_from_consensus(consensus_seq_list, junc_pos_list, junction)

        else:
            continue

        # save peptides for each isoform for each junction
        peptide_df = pd.DataFrame({'isoform': frame_list,
                                   'isoform_num': range(len(consensus_seq_list)),
                                   'junc_aa_pos': junc_aa_list, 
                                   'peptides': peptides, 
                                   'junction': junction,
                                   'frame_selection': frame_selection})
        if len(total_peptide_df)==0:
            total_peptide_df = peptide_df.copy()
        else:
            total_peptide_df = total_peptide_df.append(peptide_df)

    savepath = os.path.join(args.intermediate_results, 'neopeptides.tsv')
    total_peptide_df.to_csv(savepath, sep='\t', index=False)

