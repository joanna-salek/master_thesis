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


def distance_tree(matrix, organisms, title):
    x = pdist(matrix)
    dist = linkage(x, method="average")  # average = UPGMA method
    plt.figure(figsize=(20, 5))
    dendrogram(dist, labels=organisms, orientation="left")
    plt.title(title)
    plt.show()
