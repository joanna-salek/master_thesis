from classic_free_aligment.kmers.kmers import fasta_tree_kmers


################################################
# EXAMPLE 1
################################################
def example1():
    fasta_tree_kmers("../../in/influenza_viruses",
                     "Influenza virus phylogenetic tree - kmer method")


################################################
# EXAMPLE 2
################################################
def example2():
    fasta_tree_kmers("../../in/SPARC_refseq_transcript.fasta",
                     "SPARC gene phylogenetic tree - kmer method")
