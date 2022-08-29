################################################
# EXAMPLE MIXED
################################################
from ULTIMATE_EXAMPLE import hi, mom_inte, di
from app.helpers import get_species_list_SPARC


def MIXED(file, chart_name, organisms):
    hurst = hi(file, chart_name, organisms)
    interia = mom_inte(file, chart_name, organisms)
    DFT = di(file, chart_name, organisms)
    for r in range(len(hurst)):
        for ele in interia[r]:
            DFT[r].append(ele)
            DFT[r].append(hurst[r])
    matrix_and_tree(DFT, organisms, chart_name)


## FILES PATHS
i_file = r"C:\Users\joann\PycharmProjects\magisterka\in\influenza_viruses"
SPARC_file = r"C:\Users\joann\PycharmProjects\magisterka\in\SPARC_refseq_transcript.fasta"
C_file = r"C:\Users\joann\PycharmProjects\magisterka\in\CORONAVIRUSES1.txt"

## CHART NAMES
i_name = "Drzewo filogenetyczne wirusów grypy - gen neuroamidazy (UPGMA) - "
SPARC_name = "Drzewo filogenetyczne ssaków - gen SPARC (UPGMA) - "
c_name = "Drzewo filogenetyczne koronawirusów - cały genom (UPGA) - "

## METHOD NAMES
mixed_name = "Metoda mieszana"



MIXED(SPARC_file, SPARC_name + mixed_name, get_species_list_SPARC())
