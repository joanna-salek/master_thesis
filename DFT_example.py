from chaos_game import DFT_CGR, draw_tree, DFT_from_fasta, T_m
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

      d = DFT_CGR(seq)
      print(d.get_DFT())
      d.plot_CGR()

################################################
# EXAMPLE 2
################################################
# read file from fasta with only one seq in it
# you can enable drawing chaos game representation
# and choose CGR type
def example2():
      dft = DFT_from_fasta("one_file_test.fasta", True, CGR_types="RY")
      for key in dft.keys():
            print(key, dft[key])


################################################
# EXAMPLE 3
################################################
# read file from fasta with multiple seq in it
# fasta parser will divide it on sequences and species
def example3():
      dft = DFT_from_fasta("SPARC_refseq_transcript.fasta")
      for key in dft.keys():
            print(key, dft[key])


################################################
# EXAMPLE 4
################################################
def example4():
      dft = DFT_from_fasta("SPARC_refseq_transcript.fasta")
      arr = []
      for key in dft["WS"].keys():
            arr.append(dft["WS"][key])

      n = max([len(x) for x in arr])  # maksymalna dlugosc sewkecni
      arr = [T_m(x, n)[1:] for x in arr]  # wydluzone widma (do najdluzszej)

      draw_tree(arr, list(dft["WS"].keys()), "SPARC gene phylogenetic tree - DFT method")

# example1()
# example2()
# example3()
# example4()