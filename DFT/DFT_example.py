from DFT import DFT_CGR, DFT_from_fasta, T_m, DFT_data
from app.helpers import draw_tree, euclidean_matrix, distance_tree, get_species_list_coronaviruses, \
    get_NCBI_IDS_list_influenza, get_species_list_SPARC

################################################
# EXAMPLE 1
################################################
# read seq directly by pasting it
from app.read_data import fasta_parser


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
    dft = DFT_from_fasta("../in/one_file_test.fasta", True, CGR_types="RY")
    for key in dft.keys():
        print(key, dft[key])


################################################
# EXAMPLE 3
################################################
# read file from fasta with multiple seq in it
# fasta parser will divide it on sequences and species
def example3():
    dft = DFT_from_fasta("../in/SPARC_refseq_transcript.fasta")
    for key in dft.keys():
        print(key, dft[key])


################################################
# EXAMPLE 4
################################################
def example4():
    dft = DFT_from_fasta("../in/SPARC_refseq_transcript.fasta")
    print(dft)
    arr = []
    for key in dft["WS"].keys():
        arr.append(dft["WS"][key])
    n = max([len(x) for x in arr])  # maksymalna dlugosc sewkecni
    arr = [T_m(x, n)[1:] for x in arr]  # wydluzone widma (do najdluzszej)
    draw_tree(arr, list(dft["WS"].keys()), "SPARC gene phylogenetic tree - DFT method")


################################################
# EXAMPLE 5
################################################
def example5():
    data = fasta_parser("../in/influenza_viruses")
    arr = []
    for d in data[1]:
        arr.append(DFT_data(d))
    n = max([len(x) for x in arr])  # maksymalna dlugosc sewkecni
    arr = [T_m(x, n)[1:] for x in arr]  # wydluzone widma (do najdluzszej)
    matrix = euclidean_matrix(arr)
    print (get_NCBI_IDS_list_influenza())
    print (len(matrix))

    distance_tree(matrix, get_NCBI_IDS_list_influenza(), "Drzewo filogenetyczne wirusów grypy - Metoda DFT (Odległość euklidesowa)")


################################################
# EXAMPLE 6
################################################
def example6():
    data = fasta_parser("../in/coronaviruses.txt")
    arr = []
    for d in data[1]:
        arr.append(DFT_data(d))
    n = max([len(x) for x in arr])  # maksymalna dlugosc sewkecni
    arr = [T_m(x, n)[1:] for x in arr]  # wydluzone widma (do najdluzszej)
    matrix = euclidean_matrix(arr)
    title = "Drzewo filogenetyczne Koronawirusów - Metoda DFT (Odległość euklidesowa)"
    distance_tree(matrix, get_species_list_coronaviruses(), title)

################################################
# EXAMPLE 7
################################################
def example7():
    data = fasta_parser("../in/SPARC_refseq_transcript.fasta")
    arr = []
    for d in data[1]:
        arr.append(DFT_data(d))
    n = max([len(x) for x in arr])  # maksymalna dlugosc sewkecni
    arr = [T_m(x, n)[1:] for x in arr]  # wydluzone widma (do najdluzszej)
    matrix = euclidean_matrix(arr)
    title = "Drzewo filogenetyczne SPARC - Metoda DFT (Odległość euklidesowa)"
    distance_tree(matrix, get_species_list_SPARC(), title)


example7()