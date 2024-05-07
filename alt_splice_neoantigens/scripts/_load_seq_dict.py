from collections import defaultdict
from Bio import SeqIO

# load CDS/cDNA file

def load_CDS_dict(CDS_path):
    CDS_dict = defaultdict(list)
    strand_dict = defaultdict(list)
    gene_dict = defaultdict(list)
    with open(CDS_path, 'r') as handle:
        for seq_record in SeqIO.parse(handle, 'fasta'):
            enst = seq_record.id.split('.')[0]
            if enst in CDS_dict:
                # check which one is more accurate
                print(enst) #### none are duplicates!
            if str(seq_record.seq) == 'Sequenceunavailable':
                continue
            CDS_dict[enst] = str(seq_record.seq)  
            header = seq_record.description
            header = header.split(' ')
            for val in header:
                if val.startswith('chromosome'):
                    chrome = val
                    info = chrome.split(':')
                    strand_direction = info[-1]
                if val.startswith('gene_symbol'):
                    gene_sym = val
                    gene_info = gene_sym.split(':')
                    gene_name = gene_info[-1]
            gene_dict[enst] = gene_name
            strand_dict[enst] = str(strand_direction)
    return CDS_dict, strand_dict, gene_dict

