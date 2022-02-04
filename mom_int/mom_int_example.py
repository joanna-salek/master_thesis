from app.helpers import draw_tree
from mom_int import DNA, intertia_from_fasta


################################################
# EXAMPLE 1
################################################
# read seq directly by pasting it
def example1():
    seq = "AATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCT" \
          "CTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGA" \
          "GAGAGAAATCTTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTC" \
          "TGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAG" \
          "AGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAATCTTCTCT" \
          "GAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGA" \
          "GAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATC" \
          "TTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAG" \
          "AGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCT" \
          "TCTCTGAGAGAGAGAGAAATCTTCTCTGAGAG"

    i = DNA(seq, 2, False)
    print(i.DNA_descriptors())


################################################
# EXAMPLE 2
################################################
# read file from fasta with only one seq in it
def example2():
    intertia = intertia_from_fasta("../in/one_file_test.fasta")
    for key in intertia.keys():
        print(key, intertia[key])


################################################
# EXAMPLE 3
################################################
# read file from fasta with multiple seq in it
# fasta parser will divide it on sequences and species
def example3():
    intertia = intertia_from_fasta("../in/SPARC_refseq_transcript.fasta")
    for key in intertia.keys():
        print(key, intertia[key])


################################################
# EXAMPLE 6
################################################
# draw phylogenetic tree based on hurst_exponent values for species
def example6():
    intertia = intertia_from_fasta("../in/SPARC_refseq_transcript.fasta")
    arr = []
    for key in intertia.keys():
        arr.append(intertia[key])
    draw_tree(arr, list(intertia.keys()), "SPARC gene phylogenetic tree - moments of inertia")


################################################
# EXAMPLE 5
################################################
# Neuroamidaze influenza viruses
def example7():
    intertia = intertia_from_fasta("../in/influenza_viruses")
    arr = []
    for key in intertia.keys():
        arr.append(intertia[key])
    draw_tree(arr, list(intertia.keys()), "Influenza viruses phylogenetic tree - moments of inertia")

# example1()
# example2()
# example3()
# example4()
example6()
# example7()
