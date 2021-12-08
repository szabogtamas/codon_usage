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
from Bio import Entrez, SeqIO
Entrez.email = ""
```

```python
# Use T7 phage genome as test, inspired by https://wilkelab.org/classes/SDS348/2020_spring/worksheets/class20_solutions.html
genome_id = "NC_001604"
```

```python
# Download entire genome:
download_handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb", retmode="text")
data = download_handle.read()
download_handle.close()

out_handle = open(genome_id+".gb", "w")
out_handle.write(data)
out_handle.close()
```

```python
in_handle = open(genome_id+".gb", "r")
record = SeqIO.read(in_handle, "genbank")
in_handle.close()

for feature in record.features:
    if feature.type == 'CDS':
        locus_tag = feature.qualifiers['locus_tag'][0]
        print(feature)
        DNA_seq = target_feature.extract(record).seq
        print("DNA:", str(DNA_seq))
        print("protein:", DNA_seq.translate())
```

```python

```
