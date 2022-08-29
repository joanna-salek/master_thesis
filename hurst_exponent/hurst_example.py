from app.helpers import draw_tree, euclidean_matrix, distance_tree, get_species_list_SPARC, get_NCBI_IDS_list_influenza, \
      get_species_list_coronaviruses
from hurst_exponent.h import Hurst_CGR
from hurst_exponent import h

from app.read_data import fasta_parser


################################################
# EXAMPLE 1
################################################
# read seq directly by pasting it


def example1():
      seq = "AATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCT" \
            "CTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGA" \
            "GAGAGAAATCTTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTC" \
            "TGAGAGAGAGAGAAATCGTGAGTGCGCTAGCGATGCCGAGACCGATAGCGAT" \
            "TTCTCTGAGAGAGCGTAGGCTAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCT" \
            "CTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGA" \
            "GAGAGAAATCTTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTC" \
            "TGAGAGAGAGAGAAATCGTGAGTGCGCTAGCGATGCCGAGACCGATAGCGAT" \
            "TTCTCTGAGAGAGCGTAGGCTCGAACGCTGACGAAGATCGCTATAGCATGCG" \
            "AGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAATCTTCTCT" \
            "GAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGA" \
            "GAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATC" \
            "TTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAG" \
            "AGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCT" \
            "TCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGCGAACGCTGACGAAGATCGCTATAGCATGCG" \
            "AGAGAAATCTTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCT" \
            "CTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGA" \
            "GAGAGAAATCTTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTC" \
            "TGAGAGAGAGAGAAATCGTGAGTGCGCTAGCGATGCCGAGACCGATAGCGAT" \
            "TTCTCTGAGAGAGCGTAGGCTCGAACGCTGACGAAGATCGCTATAGCATGCG" \
            "AGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAATCTTCTCT" \
            "GAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGA" \
            "GAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATC" \
            "TTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAG" \
            "AGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCT" \
            "TCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAATCTTCTCT" \
            "GAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGA" \
            "GAGAAATCTTCTCTGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCT" \
            "CTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGA" \
            "GAGAGAAATCTTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTC" \
            "TGAGAGAGAGAGAAATCGTGAGTGCGCTAGCGATGCCGAGACCGATAGCGAT" \
            "TTCTCTGAGAGAGCGTAGGCTCGAACGCTGACGAAGATCGCTATAGCATGCG" \
            "AGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAATCTTCTCT" \
            "GAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGA" \
            "GAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATC" \
            "TTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAG" \
            "AGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCT" \
            "TCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGGAGAGAAATCTTCTCTGAGAGAGAGAGAAATC" \
            "TTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAG" \
            "AGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCT" \
            "TCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCT" \
            "CTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGA" \
            "GAGAGAAATCTTCTAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCT" \
            "CTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGA" \
            "GAGAGAAATCTTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTC" \
            "TGAGAGAGAGAGAAATCGTGAGTGCGCTAGCGATGCCGAGACCGATAGCGAT" \
            "TTCTCTGAGAGAGCGTAGGCTCGAACGCTGACGAAGATCGCTATAGCATGCG" \
            "AGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAATCTTCTCT" \
            "GAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGA" \
            "GAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATC" \
            "TTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAG" \
            "AGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCT" \
            "TCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTC" \
            "TGAGAGAGAGAGAAATCGTGAGTGCGCTAGCGATGCCGAGACCGATAGCGAT" \
            "TTCTCTGAGAGAGCGTAGGCTCGAACGCTGACGAAGATCGCTATAGCATGCG" \
            "AGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAATCTTCTCT" \
            "GAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGA" \
            "GAGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATC" \
            "TTCTCTGAGAGAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAG" \
            "AGAAATCTTCTCTGAGAGAGAGAGAAATCTTCTCTGAGAGAGAGAGAAATCT" \
            "TCTCTGAGAGAGAGAGAAATCTTCTCTGAGAG"
      z = Hurst_CGR(seq)
      print(z.get_hurst())
      z.plot_CGR()

example1()
################################################
# EXAMPLE 2
################################################
# read file from fasta with only one seq in it
# you can enable drawing chaos game representation
# and choose CGR type
def example2():
      hurst = h.hurst_from_fasta("../in/one_file_test.fasta", True)
      for key in hurst.keys():
            print(key, hurst[key])

################################################
# EXAMPLE 3
################################################
# read file from fasta with multiple seq in it
# fasta parser will divide it on sequences and species
def example3():
      hurst = h.hurst_from_fasta("../in/SPARC_refseq_transcript.fasta")
      for key in hurst.keys():
            print(key, hurst[key])

################################################
# EXAMPLE 4
################################################
# save as tables from pandas to hmtl file
def example4():
      hurst = h.hurst_from_fasta("../in/SPARC_refseq_transcript.fasta")
      h.save_hurst_table(hurst)

################################################
# EXAMPLE 5
################################################
# return dictionary
def example5():
      hurst = h.hurst_from_fasta("../in/SPARC_refseq_transcript.fasta")
      for key in hurst.keys():
            i = hurst[key]
            new_dict = {k: v for k, v in sorted(i.items(), key=lambda item: item[1])}
      return new_dict

################################################
# EXAMPLE 6
################################################
# draw phylogenetic tree based on hurst_exponent values for species
def example6():
      hurst = h.hurst_from_fasta("../in/SPARC_refseq_transcript.fasta")
      arr = []
      for key in hurst["RY"].keys():
            arr.append([hurst["RY"][key]])
      draw_tree(arr, list(hurst["RY"].keys()), "SPARC gene phylogenetic tree - Hurst")

################################################
# EXAMPLE 7
################################################
# neuroamidaza wirusów grypy
def example7():
      from statistics import mean, stdev
      hurst = h.hurst_from_fasta("../in/influenza_viruses")
      type_cgr = "WS"
      arr = []
      count_H1N1 = []
      count_H7N3 = []
      count_H5N1 = []
      count_H2N2 = []
      count_H7N9 = []
      for key in hurst[type_cgr].keys():
            arr.append([hurst[type_cgr][key]])
            type_virus = key[-5:-1]
            if type_virus == "H1N1":
                  count_H1N1.append(hurst[type_cgr][key])
            elif type_virus == "H7N3":
                  count_H7N3.append(hurst[type_cgr][key])
            elif type_virus == "H5N1":
                  count_H5N1.append(hurst[type_cgr][key])
            elif type_virus == "H2N2":
                  count_H2N2.append(hurst[type_cgr][key])
            elif type_virus == "H7N9":
                  count_H7N9.append(hurst[type_cgr][key])
      file = open(f"Influenza_hurst_statistics_{type_cgr}.txt", "a")
      file.write("H1N1 mean: " + str(mean(count_H1N1)) + " stdev: " +str(stdev(count_H1N1)) + "\n")
      file.write("H7N3 mean: " + str(mean(count_H7N3)) + " stdev: " +str(stdev(count_H7N3)) + "\n")
      file.write("H1N1 mean: " + str(mean(count_H5N1)) + " stdev: " +str(stdev(count_H5N1)) + "\n")
      file.write("H1N1 mean: " + str(mean(count_H2N2)) + " stdev: " +str(stdev(count_H2N2)) + "\n")
      file.write("H1N1 mean: " + str(mean(count_H7N9)) + " stdev: " +str(stdev(count_H7N9)) + "\n")
      file.close()
      draw_tree(arr, list(hurst[type_cgr].keys()), "Influenza phylogenetic tree - hurst_exponent method")


################################################
# EXAMPLE 8
################################################
# macież różnic odległości eukolidesowej hurst_exponent
def example8():
      data = fasta_parser("../in/influenza_viruses")
      hurst_list = []
      print (data)
      for h in data[1]:
            hurst = Hurst_CGR(h).get_hurst()
            hurst_list.append(hurst)
      matrix = euclidean_matrix(hurst_list)
      title = "influenza viruses phylogenetic tree - Hurst method (EUCLIDIAN DISTANCE)"
      organisms = data[0]
      distance_tree(matrix, organisms, title)

def example9():
      data = fasta_parser("../in/influenza_viruses")
      hurst_list = []
      for h in data[1]:
            hurst = Hurst_CGR(h).get_hurst()
            hurst_list.append(hurst)
      matrix = euclidean_matrix(hurst_list)
      title = "Drzewo filogenetyczne Wirusów Grypy - Metoda wykładnika Hursta (EUCLIDIAN DISTANCE)"
      distance_tree(matrix, get_NCBI_IDS_list_influenza(), title)

def example10():
      data = fasta_parser("../in/coronaviruses.txt")
      hurst_list = []
      for h in data[1]:
            hurst = Hurst_CGR(h).get_hurst()
            hurst_list.append(hurst)
      matrix = euclidean_matrix(hurst_list)
      title = "Drzewo filogenetyczne Koronawirusów - Metoda wykładnika Hursta (EUCLIDIAN DISTANCE)"
      distance_tree(matrix, get_species_list_coronaviruses(), title)

def example11():
      data = fasta_parser("../in/SPARC_refseq_transcript.fasta")
      hurst_list = []
      for h in data[1]:
            hurst = Hurst_CGR(h).get_hurst()
            hurst_list.append(hurst)
      matrix = euclidean_matrix(hurst_list)
      title = "Drzewo filogenetyczne SPARC - Metoda wykładnika Hursta (EUCLIDIAN DISTANCE)"
      distance_tree(matrix, get_species_list_SPARC(), title)

example1()
# example2()
# example3()
# example4()
# print (example5())
# example6()
# example9()

