# CoEvolution

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> CoEvolution

The goal of this project is to make a reproducible pipeline to generate Evolutionary Rate Correlations (ERCs) for any set of organisms.
Unlike prior work, our pipeline attempts to annotate all genes from a focal species, rather than focusing on a conservative set of orthologs among all of the species. Additionally, we allow gene trees to have varying topologies, which vastly decreases power to detect co-evolutionary interactions, but avoids potential issues from forcing genes into a specific tree topology. 

Most of the code used here is written in `Julia`, although several downstream analyses rely on `R` and `python`, and running extrenal tools uses shell script syntax.

This pipeline relies on several external tools, including:

1) Mafft
2) Diamond
3) iqTree2
4) seqkit

If you wish to use our shell script to download and prep protein data as well, you will also need:

1) eutils
2) ncbi genome datasets tools

It has been developed and tested in Unix, so adjustments may need to be made if you are working in MacOS. If you are using windows - please run under WSL2. If that does not make sense to you - this pipeline has sadly not been developed for you (yet - but we hope to make it as accessible as possible in the future).


To locally reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and will need to be downloaded independently (see data download section in `Pipeline.qmd`)

1. Go to `scripts/Pipeline.qmd` and follow the guide. In final version should also be output as html hosted here, but not bothering to yet.