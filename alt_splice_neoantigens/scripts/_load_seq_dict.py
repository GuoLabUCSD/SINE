from collections import defaultdict
from Bio import SeqIO

# load CDS/cDNA file

def load_CDS_dict(CDS_path):
    CDS_dict = defaultdict(list)
    with open(CDS_path, 'r') as handle:
        for seq_record in SeqIO.parse(handle, 'fasta'):
            enst = seq_record.id.split('.')[0]
            if enst in CDS_dict:
                # check which one is more accurate
                print(enst) #### none are duplicates!
            if str(seq_record.seq) == 'Sequenceunavailable':
                continue
            CDS_dict[enst] = str(seq_record.seq)  
    return CDS_dict

