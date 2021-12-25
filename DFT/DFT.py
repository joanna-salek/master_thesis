from app.chaos_game import CGR
from scipy.fftpack import fft
import cmath
import math
from app.read_data import fasta_parser


class DFT_CGR(CGR):
    def __init__(self, seq, kind="RY"):
        super().__init__(seq, kind)
        self.F = fft(self.complex)
        self.PS = []
        self.AS = []
        for x in self.F:
            self.PS.append((abs(x)) ** 2)
            self.AS.append(cmath.phase(x))

    def get_DFT(self):
        return self.PS


def T_m(seq, m):
    n = len(seq)
    if n == m:
        return seq
    seq2 = [seq[0]]
    for k in range(1, m):
        q = k * n / m
        r = math.floor(q)

        if r == 0:
            r = 1

        if q.is_integer():
            seq2.append(seq[int(q)])
        else:
            if r < n - 1:
                seq2.append(seq[r] + (q - r) * (seq[r + 1] - seq[r]))
            else:
                seq2.append(seq[r])
    return seq2


def DFT_data(kind, seq, draw_CGR):
    cgr = DFT_CGR(seq, kind)
    if draw_CGR:
        cgr.plot_CGR()
    DFT = cgr.get_DFT()
    return DFT


def DFT_from_fasta(input_file, draw_CGR=False, CGR_types=("RY", "MK", "WS")):
    items = {}
    data = fasta_parser(input_file)
    if isinstance(CGR_types, list) and len(CGR_types) > 1:
        for ele in CGR_types:
            d = {}
            for r in range(len(data[0])):
                d[data[0][r]] = DFT_data(ele, data[1][r], draw_CGR)
            items[ele] = d
    else:
        d = {}
        for r in range(len(data[0])):
            d[data[0][r]] = DFT_data(str(CGR_types), data[1][r], draw_CGR)
        items[str(CGR_types)] = d
    return items
