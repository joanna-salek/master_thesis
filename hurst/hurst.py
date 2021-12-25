from app.chaos_game import CGR
import matplotlib.pyplot as plt
from hurst import compute_Hc
import pandas as pd
import os
from app.read_data import fasta_parser


class Hurst_CGR(CGR):
    def __init__(self, seq, kind="RY"):
        super().__init__(seq, kind)
        self.hurst, self.c, self.data = compute_Hc(self.z_values, kind="change", simplified=True)

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


def hurst_data(kind, data, draw_CGR):
    cgr = Hurst_CGR(data, kind)
    if draw_CGR:
        cgr.plot_CGR()
    hurst = cgr.get_hurst()
    return hurst


def hurst_from_fasta(input_file, draw_CGR=False, CGR_types=("RY", "MK", "WS")):
    items = {}
    data = fasta_parser(input_file)
    if isinstance(CGR_types, list) and len(CGR_types) > 1:
        for ele in CGR_types:
            d = {}
            for r in range(len(data[0])):
                d[data[0][r]] = hurst_data(ele, data[1][r], draw_CGR)
            items[ele] = d
    else:
        d = {}
        for r in range(len(data[0])):
            d[data[0][r]] = hurst_data(str(CGR_types), data[1][r], draw_CGR)
        items[str(CGR_types)] = d
    return items


def save_hurst_table(hurst):
    out = ""
    for key in hurst.keys():
        item = hurst[key]
        df = pd.DataFrame.from_dict(item, orient='index', columns=["hurst value"])
        print(df)
        out = rf"out\index_{key}.html"
        with open(out, "w") as f:
            f.write(df.to_html())
    print(fr"output file directory {os.path.abspath(os.getcwd())}\{out}")
