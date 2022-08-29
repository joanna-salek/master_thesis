from app.helpers import draw_tree, get_species_list_coronaviruses, get_species_list_SPARC, get_NCBI_IDS_list_influenza
from app.read_data import fasta_parser
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



################################################
# EXAMPLE 8
################################################
def example8():
    data = fasta_parser("../in/coronaviruses.txt")
    d = {}
    for r in range(len(data[0])):
        dna = DNA(data[1][r], 2, False)
        d[data[0][r]] = dna.DNA_descriptors()
    draw_tree(list(d.values()), get_species_list_coronaviruses(), "Drzewo filogenetyczne Koronawirusów - Momenty Bezwładności (Odległość euklidesowa)")


def example9():
    data = fasta_parser("../in/influenza_viruses")
    d = {}
    for r in range(len(data[0])):
        dna = DNA(data[1][r], 2, False)
        d[data[0][r]] = dna.DNA_descriptors()
    draw_tree(list(d.values()), get_NCBI_IDS_list_influenza(), "Drzewo filogenetyczne Wirusów Grypy - Momenty Bezwładności (Odległość euklidesowa)")



def example10():
    data = fasta_parser("../in/SPARC_refseq_transcript.fasta")
    d = {}
    for r in range(len(data[0])):
        dna = DNA(data[1][r], 2, False)
        d[data[0][r]] = dna.DNA_descriptors()
    draw_tree(list(d.values()), get_species_list_SPARC(), "Drzewo filogenetyczne SPARC - Momenty Bezwładności (Odległość euklidesowa)")

# example1()
# example2()
# example3()
# example4()
# example6()
# example7
example10()

