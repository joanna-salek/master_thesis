from Bio import SeqIO
import re

def fasta_parser(fasta_sequences):
  species = []
  seq = []
  pattern = r'\b[A-Z][a-z]+\b \b[a-z]+\b'
  handle = open(fasta_sequences)
  fasta_sequences = SeqIO.parse(handle, 'fasta')
  for fasta in fasta_sequences:
    t = re.search(pattern, fasta.description)
    species.append(t.group())
    seq.append(str(fasta.seq))
  handle.close()
  return species, seq

