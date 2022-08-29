from hurst import compute_Hc

from mom_int.mom_int import intertia_from_fasta


values = intertia_from_fasta(r"C:\Users\joann\PycharmProjects\magisterka\in\SPARC_refseq_transcript.fasta")
data_to_hurst = []
for value in values.values():
    data_to_hurst.append(value[0] + value[1])


hurst = compute_Hc(data_to_hurst, kind="change", simplified=False)
print (hurst)