from classic_free_aligment.complexity_z.complexity_z import z_complex_matrix
from classic_free_aligment.kmers.kmers import mer_matrix
from hurst_exponent.h import Hurst_CGR
from app.read_data import fasta_parser
from DFT.DFT import T_m, DFT_data
from app.helpers import euclidean_matrix, distance_tree, get_species_list_coronaviruses, get_NCBI_IDS_list_influenza, \
    get_species_list_SPARC
from mom_int.mom_int import intertia_from_fasta


################################################
# EXAMPLE HURST
################################################


def hi(file, chart_name, organisms, chart_size):
    data = fasta_parser(file)
    hurst_list = []
    for h in data[1]:
        hurst = Hurst_CGR(h).get_hurst()
        hurst_list.append(hurst)
    matrix = euclidean_matrix(hurst_list)
    distance_tree(matrix, organisms, chart_name, chart_size)
    return hurst_list

################################################
# EXAMPLE MOM INT
################################################

def mom_inte(file, chart_name, organisms, chart_size):
    intertia = intertia_from_fasta(file)
    arr = []
    for key in intertia.keys():
        arr.append(intertia[key])
    matrix = euclidean_matrix(arr)
    distance_tree(matrix, organisms, chart_name, chart_size)
    return arr


################################################
# EXAMPLE DFT
################################################

def di(file, chart_name, organisms, chart_size):
    data = fasta_parser(file)
    arr = []
    for d in data[1]:
        arr.append(DFT_data(d))
    n = max([len(x) for x in arr])  # maksymalna dlugosc sewkecni
    arr = [T_m(x, n)[1:] for x in arr]  # wydluzone widma (do najdluzszej)
    matrix = euclidean_matrix(arr)
    distance_tree(matrix, organisms, chart_name, chart_size)
    return arr


################################################
# EXAMPLE MIXED METHOD (mom of int + Hurst)
################################################

def mixed(file, chart_name, organisms, chart_size):
    hurst = hi(file, chart_name, organisms, chart_size)
    interia = mom_inte(file, chart_name, organisms, chart_size)
    arr = []
    for r in range(len(hurst)):
        arr.append([hurst[r], interia[r][0], interia[r][1]])
    matrix = euclidean_matrix(arr)
    distance_tree(matrix, organisms, chart_name, chart_size)


################################################
# EXAMPLE CLASSIC FREE ALIGNMENT - KMERS
################################################

def k_mers(input_file, chart_name, organisms, chart_size):
    data = fasta_parser(input_file)
    matrix = mer_matrix(data[1])
    distance_tree(matrix, organisms, chart_name, chart_size)

################################################
# EXAMPLE CLASSIC FREE ALIGNMENT - LEMPEL-ZIV
################################################

def z_complexity(input_file, chart_name, organisms, chart_size):
    data = fasta_parser(input_file)
    matrix = z_complex_matrix(data[1])
    distance_tree(matrix, organisms, chart_name, chart_size)


## FILES PATHS
i_file = r"C:\Users\joann\PycharmProjects\magisterka\in\influenza_viruses"
SPARC_file = r"C:\Users\joann\PycharmProjects\magisterka\in\SPARC_refseq_transcript.fasta"
C_file = r"C:\Users\joann\PycharmProjects\magisterka\in\CORONAVIRUSES1.txt"

## CHART NAMES
i_name = "Drzewo filogenetyczne wirusów grypy - gen neuroamidazy (UPGMA) - "
SPARC_name = "Drzewo filogenetyczne ssaków - gen SPARC (UPGMA) - "
c_name = "Drzewo filogenetyczne koronawirusów - cały genom (UPGMA) - "

## METHOD NAMES
h_name = "Wykładnik Hursta"
DFT_name = "Transformata Fouriera"
mom_int_name = "Momenty bezwładności"
mixed_name = "Metoda mieszana"
classic_free_aligment = "Klasyczne metody bez-dopasowania "
lempel_z = "Zlozoność Lempel-Ziv"
k_mers_name = "K-mery"
classic_aligment = "Klasyczne metody z dopasowaniem "


# mom_inte(C_file, c_name + mom_int_name, get_species_list_coronaviruses(), (13, 12))
# m = mom_inte(i_file, i_name + mom_int_name, get_NCBI_IDS_list_influenza(), (9, 17))
# mom_inte(SPARC_file, SPARC_name + mom_int_name, get_species_list_SPARC(), (10, 12))

# hi(C_file, c_name + h_name, get_species_list_coronaviruses(), (13, 12))
# t  = hi(i_file, i_name + h_name, get_NCBI_IDS_list_influenza(), (9, 17))
# hi(SPARC_file, SPARC_name + h_name, get_species_list_SPARC(), (10, 12))

# di(C_file, c_name + DFT_name, get_species_list_coronaviruses(), (13, 12))
# di(i_file, i_name + DFT_name, get_NCBI_IDS_list_influenza(), (9, 17))
# di(SPARC_file, SPARC_name + DFT_name, get_species_list_SPARC(), (10, 12))

# z_complexity(C_file, c_name + classic_free_aligment + lempel_z, get_species_list_coronaviruses(), (16, 12))
# z_complexity(SPARC_file, SPARC_name + classic_free_aligment + lempel_z, get_species_list_SPARC(), (12, 12))
# z_complexity(i_file, i_name + lempel_z, get_NCBI_IDS_list_influenza(), (10, 17))

# k_mers(C_file, c_name + classic_free_aligment + k_mers_name, get_species_list_coronaviruses(), (14.5, 12))
# k_mers(SPARC_file, SPARC_name + classic_free_aligment + k_mers_name, get_species_list_SPARC(), (12, 12))
# k_mers(i_file, i_name + classic_free_aligment + k_mers_name, get_NCBI_IDS_list_influenza(), (10.7, 17))

# mixed(SPARC_file, SPARC_name + mixed_name, get_species_list_SPARC(), (10, 12))
# mixed(C_file, c_name + mixed_name, get_species_list_coronaviruses(), (13, 12))
# mixed(i_file, i_name + mixed_name, get_NCBI_IDS_list_influenza(), (9, 17))


