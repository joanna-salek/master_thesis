import numpy as np

from app.helpers import distance_tree
from app.read_data import fasta_parser


def word_seq(seq):
    data = []
    length = len(seq)
    i = 0
    k = 1
    while i < length:
        while seq[i:i + k] in data and i + k < length:
            k += 1
        if seq[i:i + k] not in data:
            data.append(seq[i:i + k])
        i += k
        k = 1
    return data


def z_complexity(seq1, seq2):
    N1 = len(word_seq(seq1))
    N2 = len(word_seq(seq2))
    N = len(word_seq(seq1 + seq2))
    C = (N - min(N1, N2)) / max(N1, N2)
    return C


def z_complex_matrix(sequences):
    n = len(sequences)
    a = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            a[i][j] = z_complexity(sequences[i], sequences[j])
    return a


def fasta_tree_z_complexity(input_file, title):
    data = fasta_parser(input_file)
    matrix = z_complex_matrix(data[1])
    organisms = data[0]
    distance_tree(matrix, organisms, title)
