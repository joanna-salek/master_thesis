import unittest
import pandas as pd
from CGR import CGR
from read_data import fasta_parser

class CGR_tests(unittest.TestCase):

    def test_hurst(self):
        data = fasta_parser("SPARC_refseq_transcript.fasta")
        species = [data[0]]
        hurst = []
        CGR_types = ["RY", "MK", "WS"]
        for ele in CGR_types:
            hurst.append(create_data(ele, data[1]))
        df = pd.DataFrame(hurst, columns=species, index=CGR_types)
        with open("tables/SPARC_Hurst.html", "w") as f:
            f.write(df.to_html())
        print(df)

def create_data(kind, data):
    hurst = []
    for value in data:
        cgr = CGR(value, kind)
        hurst.append(cgr.get_hurst())
    return hurst

if __name__ == '__main__':
    unittest.main()