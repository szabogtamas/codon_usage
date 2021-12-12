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

import os
```

```python
genome_dir = "../fasta_transcriptomes"
```

```python
def create_codonfeq_table(transcriptome):
    for seq_record in SeqIO.parse(transcriptome, "fasta"):
        print(seq_record.id)
        print(seq_record.translate())
        print(repr(seq_record.seq))
        break

cft = create_codonfeq_table(genome_dir + "/" + os.listdir(genome_dir)[0])
```

```python

```

```python

```
