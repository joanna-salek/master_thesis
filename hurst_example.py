from chaos_game import Hurst_CGR, hurst_from_fasta, draw_tree, save_hurst_table

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

      h = Hurst_CGR(seq)
      print(h.get_hurst())
      h.plot_CGR()

################################################
# EXAMPLE 2
################################################
# read file from fasta with only one seq in it
# you can enable drawing chaos game representation
# and choose CGR type
def example2():
      hurst = hurst_from_fasta("one_file_test.fasta", True, CGR_types="RY")
      for key in hurst.keys():
            print(key, hurst[key])

################################################
# EXAMPLE 3
################################################
# read file from fasta with multiple seq in it
# fasta parser will divide it on sequences and species
def example3():
      hurst = hurst_from_fasta("SPARC_refseq_transcript.fasta")
      for key in hurst.keys():
            print(key, hurst[key])

################################################
# EXAMPLE 4
################################################
# save as tables from pandas to hmtl file
def example4():
      hurst = hurst_from_fasta("SPARC_refseq_transcript.fasta")
      save_hurst_table(hurst)

################################################
# EXAMPLE 5
################################################
# return dictionary
def example5():
      hurst = hurst_from_fasta("SPARC_refseq_transcript.fasta")
      for key in hurst.keys():
            i = hurst[key]
            new_dict = {k: v for k, v in sorted(i.items(), key=lambda item: item[1])}
      return new_dict

################################################
# EXAMPLE 6
################################################
# draw phylogenetic tree based on hurst values for species
def example6():
      hurst = hurst_from_fasta("SPARC_refseq_transcript.fasta")
      arr = []
      for key in hurst["MK"].keys():
            arr.append([hurst["MK"][key]])
      draw_tree(arr, list(hurst["MK"].keys()), "SPARC gene phylogenetic tree - moments of inertia")

################################################
# EXAMPLE 7
################################################
# neuroamidaza wirusów grypy
def example7():
      hurst = hurst_from_fasta("influenza_viruses")
      arr = []
      for key in hurst["WS"].keys():
            arr.append([hurst["WS"][key]])
      draw_tree(arr, list(hurst["WS"].keys()), "Influenza phylogenetic tree - hurst method")




# example1()
# example2()
# example3()
# example4()
# print (example5())
# example6()
# example7()