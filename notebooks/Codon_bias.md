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
    display_name: Python 3
    language: python
    name: python3
---

# Draft a function drawing faceted subnetworks


## Dependencies

```python
### Tools to be used
import matplotlib
from matplotlib import pyplot as plt

import numpy as np
import pandas as pd
```

```python
### Download reference sequences

!wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.align.gz
!gunzip hg38.fa.align.gz
```

```python
### Calculate codon stats

from Bio.SeqUtils import CodonUsage

usage_tab = CodonUsage.CodonAdaptationIndex
usage_tab.generate_index("hg38.fa.align")
usage_tab.print_index()
```