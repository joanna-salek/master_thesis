import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage  # to generate a tree
from scipy.spatial.distance import pdist, euclidean


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
            CGR = ([last_ele[0] - ((last_ele[0] - self.vertices[self.seq[r]][0]) * 0.5),
                    last_ele[1] - ((last_ele[1] - self.vertices[self.seq[r]][1]) * 0.5)])
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
        ax.plot(1, 0, marker='>', transform=ax.get_yaxis_transform(), **arrow_fmt)
        ax.plot(0, 1, marker='^', transform=ax.get_xaxis_transform(), **arrow_fmt)

        # Add nucleotide labels
        ax.text(1.02, 1.03, get_key(self.vertices, [1, 1]), fontsize=16, fontweight='bold', va='top')
        ax.text(-0.02, 1.03, get_key(self.vertices, [0, 1]), fontsize=16, fontweight='bold', va='top')
        ax.text(-0.03, -0.02, get_key(self.vertices, [0, 0]), fontsize=16, fontweight='bold', va='top')
        ax.text(1.02, -0.02, get_key(self.vertices, [1, 0]), fontsize=16, fontweight='bold', va='top')

        plt.show()
