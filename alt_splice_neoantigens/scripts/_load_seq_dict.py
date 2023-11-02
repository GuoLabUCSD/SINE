from collections import defaultdict
from Bio import SeqIO

# load CDS/cDNA file
#def load_CDS_dict(CDS_path='/cellar/users/andreabc/Data/refs/Homo_sapiens.GRCh38.cds.all.fa'):
#    CDS_dict = defaultdict(list)
#    with open(CDS_path, 'r') as handle:
#        for seq_record in SeqIO.parse(handle, 'fasta'):
#            enst = seq_record.id.split('.')[0]
#            if enst in CDS_dict:
#                # check which one is more accurate
#                print(enst) #### none are duplicates!
#            if str(seq_record.seq) == 'Sequenceunavailable':
#                continue
#            CDS_dict[enst] = str(seq_record.seq)  
#    return CDS_dict

def load_mouse_CDS_dict(CDS_path):
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

#def load_peptide_dict(path='/cellar/users/andreabc/Data/refs/Homo_sapiens.GRCh38.pep.all.fa', keytype='ENST'):
#    peptide_dict = defaultdict(list)
#    with open(path, 'r') as handle:
#        for seq_record in SeqIO.parse(handle, 'fasta'):
#            if keytype=='ENST':
#                key = seq_record.description.split('transcript:')[1].split('.')[0]
#            elif keytype=='gene_symbol':
#                key = seq_record.description.split('gene_symbol:')[1].split(' ')[0] 
#            if key in peptide_dict:
#                # take longest one
#                if len(peptide_dict[key])> len(str(seq_record.seq)):
#                    continue
#            if str(seq_record.seq) == 'Sequenceunavailable':
#                continue
#            peptide_dict[key] = str(seq_record.seq)  
#    return peptide_dict

