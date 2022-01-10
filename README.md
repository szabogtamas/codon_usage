# codon_usage

Explore codon frequencies in individual genes of a given reference genome.

## Quickstart

To check out pathways that are affected by unbalanced usage of a given codon,
check out the [R notebook example](notebooks/Bias_score.Rmd) that uses [precalculated
codon frequencies](stock_data/codon_freqs) available in this repo.  
You can open this in a docker container provided by the Rocker project, with
preinstalled packages, running the following command, optimally in the repo folder:  

```
docker run --rm   -p 127.0.0.1:8989:8989 -p 127.0.0.1:8787:8787 -e USERID=$UID -e PASSWORD=SecurePassword -v $PWD:/home/rstudio/local_files szabogtamas/codon_usage
```

An RStudio instance will be available at http://localhost:8787/ for your convenience.
The user name is rstudio and the password is what you typed instead of "SecurePassword"
above.

Some example figures can be seen in the [example_outputs](example_outputs) folder.

Project is under development. Some alternative approaches are being tested in the
dev branch.