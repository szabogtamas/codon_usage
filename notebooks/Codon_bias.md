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

# Create codon usage table from a reference

## Dependencies

```python
### Tools to be used
import matplotlib
from matplotlib import pyplot as plt

import numpy as np
import pandas as pd
```

```python
from Bio import SeqIO
```

```python
### Download reference sequences

#https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
#https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz
!wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.align.gz
!gunzip hg38.fa.align.gz
```

```python
for seq_record in SeqIO.parse("GRCh38_latest_rna.fna", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    break
```

```python

```

```python

```

```python

```

```python
### Calculate codon stats

from Bio.SeqUtils import CodonUsage

usage_tab = CodonUsage.CodonAdaptationIndex
usage_tab.generate_index("hg38.fa.align")
usage_tab.print_index()
```

```python
### Show frequency of a given codon

CODON = "TTT"

print(usage_tab[CODON])
```
