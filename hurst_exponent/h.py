from hurst import compute_Hc
from app.chaos_game import CGR
import matplotlib.pyplot as plt
import pandas as pd
import os
from app.read_data import fasta_parser


class Hurst_CGR(CGR):
    def __init__(self, seq, kind="RY"):
        super().__init__(seq, kind)
        self.hurst, self.c, self.data = compute_Hc(self.z_values, kind="change", simplified=False)

    def get_hurst(self):
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


def hurst_data(data, draw_CGR):
    cgr = Hurst_CGR(data)
    if draw_CGR:
        cgr.plot_CGR()
    hurst = cgr.get_hurst()
    return hurst


def hurst_from_fasta(input_file, draw_CGR=False):
    data = fasta_parser(input_file)
    d = {}
    for r in range(len(data[0])):
        d[data[0][r]] = hurst_data(data[1][r], draw_CGR)
    else:
        d = {}
        for r in range(len(data[0])):
            d[data[0][r]] = hurst_data(data[1][r], draw_CGR)
    return d


def save_hurst_table(hurst):
    out = ""
    for key in hurst.keys():
        item = hurst[key]
        df = pd.DataFrame.from_dict(item, orient='index', columns=["hurst_exponent value"])
        print(df)
        out = rf"out\index_{key}.html"
        with open(out, "w") as f:
            f.write(df.to_html())
    print(fr"output file directory {os.path.abspath(os.getcwd())}\{out}")

