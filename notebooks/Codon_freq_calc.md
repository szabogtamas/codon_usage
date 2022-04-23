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
genome_dir = "../stock_data/fasta_transcriptomes"
freq_dir = "../stock_data/codon_freqs"
pep_dir = "../stock_data/pentapep_freqs"
```

```python
def codon_freq_of_cds(seq, proto_d={"".join(x): 0 for x in itertools.product("CUAG", repeat=3)}, complement=False):
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

#cft = create_codonfeq_table(genome_dir + "/" + os.listdir(genome_dir)[0])
#cft.head()
```

```python
def create_pepfreq_table(transcriptome, peps):
    frqs = dict()
    for seq_record in SeqIO.parse(transcriptome, "fasta"):
        frqs[seq_record.id] = {x: seq_record.translate().seq.count(x) for x in peps}
    return pd.DataFrame.from_dict(frqs, orient="index")

#cft = create_codonfeq_table(genome_dir + "/" + os.listdir(genome_dir)[0])
#cft.head()
```

```python
def ffc(transcriptome):
    for seq_record in SeqIO.parse(transcriptome, "fasta"):
        #print(dir(seq_record))
        #if str(seq_record.name).find("901") > -1:
        if seq_record.name == "18119_901":
            print(seq_record.name)
            print(dir(seq_record))
            print(seq_record.seq)
            print(seq_record.dbxrefs)
            print(seq_record.description)
            print(seq_record.reverse_complement())
            print(seq_record.translate().seq.count("EERRR"))
            print(seq_record.translate().seq)
    return

ffc(genome_dir + "/" + os.listdir(genome_dir)[-1])
```

```python

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
for fn in os.listdir(genome_dir):
    freq_tab = create_pepfreq_table(genome_dir + "/" + fn, ["EERRR", "EEERR", "EERRE"])
    freq_tab.to_csv(pep_dir + "/" + os.path.splitext(fn)[0] +"_pentapep.csv")
```

```python
df = pd.read_csv(pep_dir + "/" + os.path.splitext(fn)[0] +"_pentapep.csv")
df = df.loc[df.EERRR > 0,:]
df.shape
```

```python

```
