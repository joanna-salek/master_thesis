from io import StringIO
import matplotlib.pyplot as plt
from Bio import Phylo
from hurst import compute_Hc
from read_data import fasta_parser
import os
import pandas as pd
import cmath
from scipy.fftpack import fft
import math
from scipy import spatial
import numpy as np



class CGR:

    def __init__(self, seq, kind="RY"):
        self.coordinates = []
        self.z_values = []
        self.vertices = {}
        self.complex = []
        self.seq = seq
        if kind in {"RY", "MK", "WS"}:
            self.kind = kind
        else:
            self.kind = "RY"
            print("Warning: not valid kind of vertices. Vertices set to defaulft RY")
        self.assign_vertices()
        self.calculate_CGR()

    def get_z_values(self):
        return self.z_values

    def get_coordinates(self):
        return self.coordinates

    def get_complex(self):
        return self.complex

    def assign_vertices(self):
        if self.kind == "RY":
            self.vertices["A"] = [0, 0]
            self.vertices["C"] = [0, 1]
            self.vertices["G"] = [1, 1]
            self.vertices["T"] = [1, 0]
        elif self.kind == "MK":
            self.vertices["A"] = [0, 0]
            self.vertices["C"] = [1, 1]
            self.vertices["G"] = [0, 1]
            self.vertices["T"] = [1, 0]
        elif self.kind == "WS":
            self.vertices["A"] = [0, 0]
            self.vertices["C"] = [0, 1]
            self.vertices["G"] = [1, 0]
            self.vertices["T"] = [1, 1]

    def calculate_CGR(self):
        for r in range(len(self.seq)):
            if r == 0:
                last_ele = [0.5, 0.5]
            else:
                last_ele = self.coordinates[-1]
            CGR = ([last_ele[0] - ((last_ele[0] - self.vertices[self.seq[r]][0])*0.5), last_ele[1] - ((last_ele[1] - self.vertices[self.seq[r]][1])*0.5)])
            self.coordinates.append(CGR)
            self.complex.append(complex(CGR[0], CGR[1]))
            self.z_values.append(CGR[0] + CGR[1])

    def plot_CGR(self):
        xs = [ele[0] for ele in self.coordinates]
        ys = [ele[1] for ele in self.coordinates]
        color = ['g']

        # Select length of axes and the space between tick labels
        xmin, xmax, ymin, ymax = 0, 1, 0, 1

        # Plot points
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.scatter(xs, ys, c=color)

        # Set identical scales for both axes
        ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), aspect='equal')

        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Draw major and minor grid lines
        ax.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.2)

        # Draw arrows
        arrow_fmt = dict(markersize=4, color='black', clip_on=False)
        ax.plot((1), (0), marker='>', transform=ax.get_yaxis_transform(), **arrow_fmt)
        ax.plot((0), (1), marker='^', transform=ax.get_xaxis_transform(), **arrow_fmt)

        # Add nucleotide labels
        ax.text(1.02, 1.03, get_key(self.vertices, [1, 1]), fontsize=16, fontweight='bold', va='top')
        ax.text(-0.02, 1.03, get_key(self.vertices, [0, 1]), fontsize=16, fontweight='bold', va='top')
        ax.text(-0.03, -0.02, get_key(self.vertices, [0, 0]), fontsize=16, fontweight='bold', va='top')
        ax.text(1.02, -0.02, get_key(self.vertices, [1, 0]), fontsize=16, fontweight='bold', va='top')

        plt.show()


################################################
# HURST calculations
################################################
class Hurst_CGR(CGR):
    def __init__(self, seq, kind="RY"):
        super().__init__(seq, kind)

    def get_hurst(self):
        self.hurst, self.c, self.data = compute_Hc(self.z_values, kind="change", simplified=True)
        return self.hurst

    def plot_hurst(self):
        self.get_hurst()
        f, ax = plt.subplots()
        ax.plot(self.data[0], self.c * self.data[0] ** self.hurst, color="deepskyblue")
        ax.scatter(self.data[0], self.data[1], color="purple")
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Time interval')
        ax.set_ylabel('R/S ratio')
        ax.grid(True)
        plt.show()


def draw_tree(values, keys, graph_title):
    from scipy.cluster.hierarchy import dendrogram, linkage
    dist = linkage(values, method="average")
    plt.figure(figsize=(20, 5))
    plt.title(graph_title)
    dendrogram(dist, orientation="left",labels=keys, above_threshold_color='black',
                    color_threshold=150)
    plt.show()


def get_key(my_dict, val):
    for key, value in my_dict.items():
        if val == value:
            return key
    return "key doesn't exist"

def hurst_data(kind, data, draw_CGR):
    cgr = Hurst_CGR(data, kind)
    if draw_CGR == True:
        cgr.plot_CGR()
    hurst = cgr.get_hurst()
    return hurst

def hurst_from_fasta(input_file, draw_CGR=False, CGR_types=["RY", "MK", "WS"]):
    all = {}
    data = fasta_parser(input_file)
    if isinstance(CGR_types, list) and len(CGR_types) > 1:
        for ele in CGR_types:
            d = {}
            for r in range(len(data[0])):
                d[data[0][r]] = hurst_data(ele, data[1][r], draw_CGR)
            all[ele] = d
    else:
        d = {}
        for r in range(len(data[0])):
            d[data[0][r]] = hurst_data(str(CGR_types), data[1][r], draw_CGR)
        all[str(CGR_types)] = d
    return all

def save_hurst_table(hurst):
        for key in hurst.keys():
            item = hurst[key]
            df = pd.DataFrame.from_dict(item, orient='index', columns=["hurst value"])
            print (df)
            out = f"out\index_{key}.html"
            with open(out, "w") as f:
                f.write(df.to_html())
        print(fr"output file directory {os.path.abspath(os.getcwd())}\{out}")

################################################
# DFT calculations
################################################

class DFT_CGR(CGR):
    def __init__(self, seq, kind="RY"):
        super().__init__(seq, kind)

    def get_DFT(self):
        self.F = fft(self.complex)
        self.PS = []
        self.AS = []
        for x in self.F:
            self.PS.append((abs(x)) ** 2)
            self.AS.append(cmath.phase(x))
        return self.PS


def T_m(seq, m):
    n = len(seq)
    if n == m:
        return seq
    seq2 = []
    seq2.append(seq[0])
    for k in range(1, m):
        Q = k * n / m
        R = math.floor(Q)

        if R == 0:
            R = 1

        if (Q).is_integer():
            seq2.append(seq[int(Q)])
        else:
            if R < n - 1:
                seq2.append(seq[R] + (Q - R) * (seq[R + 1] - seq[R]))
            else:
                seq2.append(seq[R])
    return seq2

def DFT_data(kind, data, draw_CGR):
    cgr = DFT_CGR(data, kind)
    if draw_CGR == True:
        cgr.plot_CGR()
    hurst = cgr.get_DFT()
    return hurst

def DFT_from_fasta(input_file, draw_CGR=False, CGR_types=["RY", "MK", "WS"]):
    all = {}
    data = fasta_parser(input_file)
    if isinstance(CGR_types, list) and len(CGR_types) > 1:
        for ele in CGR_types:
            d = {}
            for r in range(len(data[0])):
                d[data[0][r]] = DFT_data(ele, data[1][r], draw_CGR)
            all[ele] = d
    else:
        d = {}
        for r in range(len(data[0])):
            d[data[0][r]] = DFT_data(str(CGR_types), data[1][r], draw_CGR)
        all[str(CGR_types)] = d
    return all


################################################
# Moments of inertia
################################################
class DNA:
    def __init__(self, seq, dim, desc1):
        self.seq = seq
        self.dim = dim
        self.N = len(seq)
        self.desc1 = desc1

    def DNA_coord(self):

        DNA_coord = np.zeros((self.N, self.dim))

        if self.dim == 2:

            DNA_dict = {'A': np.array([-1, 0]), 'C': np.array([0, 1]), 'G': np.array([1, 0]), 'T': np.array([0, -1])}

            for idx, elem in enumerate(self.seq):
                vec = DNA_dict[elem]
                DNA_coord[idx] = DNA_coord[idx - 1] + vec

        elif self.dim == 3:

            DNA_dict = {'A': np.array([-1, 0, 1]), 'C': np.array([0, 1, 1]), 'G': np.array([1, 0, 1]),
                        'T': np.array([0, -1, 1])}

            for idx, elem in enumerate(self.seq):
                vec = DNA_dict[elem]
                DNA_coord[idx] = DNA_coord[idx - 1] + vec

        return DNA_coord

    def DNA_mc(self):
        coord = self.DNA_coord()
        cm_coord = np.sum(coord, axis=0) / len(coord)

        return cm_coord

    def DNA_I_ten(self):

        center_mass = self.DNA_mc()
        coord = self.DNA_coord()
        coord_upd = np.transpose(coord - center_mass)

        I_ten = np.zeros((self.dim, self.dim))

        if self.dim == 2:
            for i in range(2):
                for j in range(2):
                    if i == 0 and j == 0:
                        I_ten[0][0] = np.sum((coord_upd[1] ** 2))
                    elif i == 1 and j == 1:
                        I_ten[1][1] = np.sum((coord_upd[0] ** 2))
                    else:
                        I_ten[i][j] = -np.sum(coord_upd[0] * coord_upd[1])

        elif self.dim == 3:
            for i in range(3):
                for j in range(3):
                    if i == 0 and j == 0:
                        I_ten[0][0] = np.sum((coord_upd[1] ** 2) + (coord_upd[2] ** 2))
                    elif i == 1 and j == 1:
                        I_ten[1][1] = np.sum((coord_upd[0] ** 2) + (coord_upd[2] ** 2))
                    elif i == 2 and j == 2:
                        I_ten[2][2] = np.sum((coord_upd[0] ** 2) + (coord_upd[1] ** 2))
                    else:
                        I_ten[i][j] = -np.sum(coord_upd[i] * coord_upd[j])

        return I_ten

    def DNA_descriptors(self):

        DNA_I = self.DNA_I_ten()
        DNA_cent_mass = self.DNA_mc()

        eig, v = np.linalg.eig(np.array(DNA_I))  # those aren't sorted correctly

        # a solution to sort the eigenvalues correctly (hope this one is a correct solution)

        d = np.sum(DNA_I, axis=0).tolist()
        w = list(enumerate(d))
        w.sort(key=lambda x: x[1])

        eig2 = [x for _, x in sorted(zip([x[0] for x in w], eig), key=lambda pair: pair[
            0])]  # the eigenvalues are now in pretty much the same order as they were in tables in the article

        # creating descriptor matrix

        D = np.zeros(self.dim)

        if self.dim == 2 and self.desc1 == False:
            for i in range(len(D)):
                D[i] = DNA_cent_mass[i] / eig2[i]

            return D

        elif self.dim == 3:
            for i in range(len(D)):
                D[i] = np.sqrt(eig2[i] / self.N)

            # Cij matrix
            C = np.zeros((3, 3))

            for i in range(3):
                for j in range(3):
                    C[i][j] = spatial.distance.cosine(DNA_cent_mass[i], eig2[j])

            return D, C


def intertia_from_fasta(input_file):
    data = fasta_parser(input_file)
    d = {}
    for r in range(len(data[0])):
        dna = DNA(data[1][r], 2, False)
        d[data[0][r]] = dna.DNA_descriptors()
    return d



