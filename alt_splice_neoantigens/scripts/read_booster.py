import pandas as pd
import sys
import os
import argparse
from Bio import SeqIO

def boost_reads(path_to_fqs, n):
    for file in os.listdir(path_to_fqs):
        if file.endswith('in.ase.fq'):
            new_l_name = file.split('.')
            left = os.path.join(path_to_fqs, file)
        elif file.endswith('in2.ase.fq'):
            new_r_name = file.split('.')
            right = os.path.join(path_to_fqs, file)
    new_l_name = os.path.join(path_to_fqs, new_l_name[0] + '.trinity_in_boosted.ase.fq')
    new_r_name = os.path.join(path_to_fqs, new_r_name[0] + '.trinity_in2_boosted.ase.fq')
    l = open(new_l_name, 'w')
    r = open(new_r_name, 'w')
    for i in range(n):
        for left_record in SeqIO.parse(left, 'fastq'):
            l.write('@' + left_record.id.split('/')[0] + '{}/1'.format('x'*i) + '\n')
            l.write(str(left_record.seq) + '\n')
            l.write('+\n')
            l.write(left_record.format('fastq').split('\n')[3] + '\n')
    l.close()
    for i in range(n):
        for right_record in SeqIO.parse(right, 'fastq'):
            r.write('@' + right_record.id.split('/')[0] + '{}/2'.format('x'*i) + '\n')
            r.write(str(right_record.seq) + '\n')
            r.write('+\n')
            r.write(right_record.format('fastq').split('\n')[3] + '\n')
    r.close()
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    path_args = parser.add_argument_group("Input options:")

    path_args.add_argument('--event_dir', type=str, help='path to event working directory to save new fastq files', required=True)
    path_args.add_argument('--read_boost', type=int, help='integer specifying the factor to duplicate reads by', required=True)

    args = parser.parse_args()

    boost_reads(args.event_dir, args.read_boost)
