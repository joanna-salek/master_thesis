import re

from Bio import SeqIO


def fasta_parser(fasta_sequences):
    species = []
    seq = []
    pattern = r'\b[A-Z][a-z]+\b \b[a-z]+\b'
    handle = open(fasta_sequences)
    fasta_sequences = SeqIO.parse(handle, 'fasta')
    for fasta in fasta_sequences:
        t = re.search(pattern, fasta.description)
        try:
            species.append(t.group())
        except AttributeError:
            species.append(fasta.description)
        seq.append((str(fasta.seq)).upper())
    handle.close()
    return species, seq
