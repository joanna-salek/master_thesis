from app.helpers import *
from app.read_data import fasta_parser


def mer_dict(seq, k=3):
    d = {}
    N = len(seq)
    for i in range(N - k + 1):
        mer = seq[i:(i + k)]
        if mer not in d.keys():
            d[mer] = 1
        else:
            d[mer] += 1
    return d


def compare_sequences(seq1, seq2, k=3):
    m1 = mer_dict(seq1, k)
    m2 = mer_dict(seq2, k)
    mers = []
    for elem in m1.keys():
        mers.append(elem)
    for elem in m2.keys():
        if elem not in mers:
            mers.append(elem)
    seq1M = []
    seq2M = []
    for elem in mers:
        if elem in m1.keys():
            seq1M.append(m1[elem])
        else:
            seq1M.append(0)

        if elem in m2.keys():
            seq2M.append(m2[elem])
        else:
            seq2M.append(0)
    d = 0
    for i in range(len(mers)):
        d += (seq1M[i] - seq2M[i]) ** 2
    d = d ** 0.5
    return d


def mer_matrix(sequences):
    n = len(sequences)
    arr = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            arr[i][j] = compare_sequences(sequences[i], sequences[j])
    return arr


def fasta_tree_kmers(input_file, title):
    data = fasta_parser(input_file)
    matrix = mer_matrix(data[1])
    organisms = data[0]
    distance_tree(matrix, organisms, title)
