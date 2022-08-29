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


def DFT_data(seq, draw_CGR=False):
    cgr = DFT_CGR(seq)
    if draw_CGR:
        cgr.plot_CGR()
    return cgr.get_DFT()


def DFT_from_fasta(input_file, draw_CGR):
    items = {}
    data = fasta_parser(input_file)
    d = {}
    for r in range(len(data[0])):
        d[data[0][r]] = DFT_data(data[1][r], draw_CGR)
    else:
        d = {}
        for r in range(len(data[0])):
            d[data[0][r]] = DFT_data(data[1][r], draw_CGR)
    return d
