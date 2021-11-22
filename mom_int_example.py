from chaos_game import DNA, intertia_from_fasta, draw_tree
import numpy as np

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
      print (i.DNA_descriptors())

################################################
# EXAMPLE 2
################################################
# read file from fasta with only one seq in it
def example2():
      intertia = intertia_from_fasta("one_file_test.fasta")
      for key in intertia.keys():
            print(key, intertia[key])

################################################
# EXAMPLE 3
################################################
# read file from fasta with multiple seq in it
# fasta parser will divide it on sequences and species
def example3():
    intertia = intertia_from_fasta("SPARC_refseq_transcript.fasta")
    for key in intertia.keys():
        print(key, intertia[key])


################################################
# EXAMPLE 5
################################################
# draw phylogenetic tree based on hurst values for species
def example6():
      intertia = intertia_from_fasta("SPARC_refseq_transcript.fasta")
      arr = []
      for key in intertia.keys():
            arr.append(intertia[key])
      draw_tree(arr, list(intertia.keys()), "SPARC gene phylogenetic tree - hurst method")

# example1()
# example2()
# example3()
# example4()
example6()