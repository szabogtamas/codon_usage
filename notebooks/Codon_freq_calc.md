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
```

```python
def create_codonfeq_table(transcriptome):
    frqs = dict()
    proto_d = {"".join(x): 0 for x in itertools.product("CUAG", repeat=3)}
    proto_d["SUM"] = 0
    for seq_record in SeqIO.parse(transcriptome, "fasta"):
        d = proto_d.copy()
        seq = seq_record.seq.reverse_complement()
        c = 0
        for i in range(0, len(seq), 3):
            c += 1
            codon = str(seq[i:i+3])
            d[codon] += 1
        d["SUM"] = c
        frqs[seq_record.id] = d
        return pd.DataFrame.from_dict(frqs, orient="index")
    return frqs

cft = create_codonfeq_table(genome_dir + "/" + os.listdir(genome_dir)[0])
cft.head()
```

```python
s = "asdfagasgdadsha"
len(s)
```

```python
s[::3]
```

```python
list(itertools.islice(s, 0, 2, 3))
```

```python
dir(itertools)
```

```python
for i in range(0, len(s), 3):
    print(s[i:i+3])
```

```python

```
