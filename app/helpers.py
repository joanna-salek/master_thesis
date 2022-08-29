import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, euclidean


def draw_tree(values, keys, graph_title):
    dist = linkage(values, method="average")
    plt.figure(figsize=(20, 5))
    plt.title(graph_title)
    dendrogram(dist, orientation="left", labels=keys, above_threshold_color='black',
               color_threshold=150)
    plt.show()


def get_key(my_dict, val):
    for key, value in my_dict.items():
        if val == value:
            return key
    return "key doesn't exist"


def euclidean_matrix(values):
    n = len(values)
    a = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            a[i][j] = euclidean(values[i], values[j])
    return a


def distance_tree(matrix, organisms, title, size):
    x = pdist(matrix)
    dist = linkage(x, method="average")  # average = UPGMA method
    plt.figure(figsize=size)
    dendrogram(dist, labels=organisms, orientation="left")
    plt.title(title)
    # plt.suptitle("Scatterplot " + str(name) + " , " + r'$\Delta$' + "Output , Zeit= " + str(time) + " s", fontsize=20, y=0.95)
    plt.show()

def get_species_list_coronaviruses():
    species = ['Alphacoronavirus/Nietoperz', 'Alphacoronavirus/Nietoperz', 'Betacoronavirus/Nietoperz', 'Betacoronavirus_MERS-CoV/Czlowiek_Rozumny', 'DeltacoronavirusSwinia_domowa', 'Deltacoronavirus/Swinia_domowa', 'Gammacoronavirus/Kura_Domowa', 'Gammacoronavirus/Ptak', 'Betacoronaviruses_MERS-CoV/Wielblad', 'Betacoronaviruses_SARS-CoV2(Alfa)/Czlowiek_Rozumny', 'Betacoronaviruses_SARS-CoV2(Delta)/Czlowiek_Rozumny', 'Betacoronaviruses_SARS-CoV2(Beta)/Czlowiek_Rozumny', 'Betacoronaviruses_SARS-CoV2(Omicron)/Czlowiek_Rozumny']
    return species

def get_NCBI_IDS_list_influenza():
    species = ["H1N1/HM370969", "H1N1/CY138562", "H1N1/CY149630", "H1N1/KC608160", "H1N1/AM157358", "H1N1/AB470663",
               "H1N1/AB546159", "H1N1/HQ897966", "H1N1/EU026046", "H1N1/FJ357114", "H1N1/GQ411894", "H1N1/CY140047",
               "H1N1/KM244078", "H5N1/HQ185381", "H5N1/HQ185383", "H5N1/EU635875", "H5N1/FM177121", "H5N1/AM914017",
               "H5N1/KF572435", "H5N1/AF509102", "H5N1/AB684161", "H5N1/EF541464", "H5N1/JF699677", "H5N1/GU186511",
               "H7N3/EU500854", "H7N3/CY129336", "H7N3/CY076231", "H7N3/CY039321", "H7N3/AY646080", "H7N9/KF259734",
               "H7N9/KF938945", "H7N9/KF259688", "H7N9/KC609801", "H7N9/CY014788", "H7N9/CY186004", "H2N2/DQ017487",
               "H2N2/CY005540", "H2N2/JX081142"]
    return species

def get_species_list_SPARC():
    species = ["Człowiek rozumny", "Mysz domowa", "Szczur wędrowny", "Bydło domowe", "Danio pręgowany",
               "Kura domowa", "Afrykańska żaba szponiasta", "Szympans zwyczajny", "Dzik euroazjatycki",
               "Makak królewski", "Królik europejski", "Ryżanka japońska", "Koń Domowy", "Wilk szary"
]
    return species


