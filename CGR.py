import matplotlib.pyplot as plt
from hurst import compute_Hc

class CGR:

    def __init__(self, seq, kind="RY"):
        self.coordinates = []
        self.z_values = []
        self.vertices = {}
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
            self.z_values.append(CGR[0] + CGR[1])


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
