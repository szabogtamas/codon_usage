from Bio import SeqIO

import numpy as np
import pandas as pd

import os, itertools

def calc_codon_freq_of_cds(
    seq,
    proto_d={"".join(x): 0 for x in itertools.product("CUAG", repeat=3)},
    complement=True
):
    c = 0
    d = proto_d.copy()
    if complement:
        seq = seq.reverse_complement()
    for i in range(0, len(seq), 3):
        c += 1
        codon = str(seq[i:i+3])
        try:
            d[codon] += 1
        except:
            pass
    d["SUM"] = c
    return d