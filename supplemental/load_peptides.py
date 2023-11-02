from collections import defaultdict
from Bio import SeqIO

def load_peptide_dict(path):
    '''Loads and parses patient.peptide file
    '''
    peptide_dict = defaultdict(list)
    with open(path, 'rU') as handle:
        for seq_record in SeqIO.parse(handle, 'fasta'):
            mut_id = '_'.join(seq_record.id.split(','))
            peptide_dict[mut_id] = str(seq_record.seq).split(',')
    
    return peptide_dict
