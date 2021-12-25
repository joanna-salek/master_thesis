from classic_free_aligment.complexity_z.complexity_z import fasta_tree_z_complexity


################################################
# EXAMPLE 1
################################################
def example1():
    fasta_tree_z_complexity("../../in/influenza_viruses.fasta",
                            "Influenza viruses phylogenetic tree - z complexity method")


################################################
# EXAMPLE 2
################################################
def example2():
    fasta_tree_z_complexity("../../in/SPARC_refseq_transcript.fasta",
                            "SPARC phylogenetic tree - z complexity method")
