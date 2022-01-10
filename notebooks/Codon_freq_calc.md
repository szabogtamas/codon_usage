---
jupyter:
  jupytext:
    formats: md,ipynb
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.0
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Calculate codon frequencies for each CDS

```python
import matplotlib
from matplotlib import pyplot as plt

from Bio import SeqIO

import numpy as np
import pandas as pd

import os, itertools
```

```python
genome_dir = "../fasta_transcriptomes"
freq_dir = "../codon_freqs"
```

```python
def codon_freq_of_cds(seq, proto_d={"".join(x): 0 for x in itertools.product("CUAG", repeat=3)}, complement=True):
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
```

```python
def create_codonfeq_table(transcriptome):
    frqs = dict()
    proto_d = {"".join(x): 0 for x in itertools.product("CUAG", repeat=3)}
    for seq_record in SeqIO.parse(transcriptome, "fasta"):
        if (len(seq_record.seq) % 3) == 0:
            frqs[seq_record.id] = codon_freq_of_cds(seq_record.seq, proto_d=proto_d)
    return pd.DataFrame.from_dict(frqs, orient="index")

cft = create_codonfeq_table(genome_dir + "/" + os.listdir(genome_dir)[0])
cft.head()
```

```python
cft.reset_index().to_csv(freq_dir + "/test_freqs.csv", index=False)
```

```python
for fn in os.listdir(genome_dir):
    freq_tab = create_codonfeq_table(genome_dir + "/" + fn)
    freq_tab.to_csv(freq_dir + "/" + os.path.splitext(fn)[0] +"_freq.csv")
```

```python

```
