# CoEvolution

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> CoEvolution

The goal of this project is to make a reproducible pipeline to generate Evolutionary Rate Correlations (ERCs) for any set of organisms.
Unlike prior work, our pipeline attempts to annotate all genes from a focal species, rather than focusing on a conservative set of orthologs among all of the species.
This pipeline relies on several external tools, including:

1) Muscle v5+
2) Diamond
3) iqTree2
4) seqkit

If you wish to use our shell script to download and prep protein data as well, you will also need:

1) eutils
2) ncbi genome datasets tools

It has been developed and tested in Unix, so adjustments may need to be made if you are working in MacOS. If you are using windows - please run under WSL2. If that does not make sense to you - this pipeline has sadly not been developed for you (yet).


To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and will need to be downloaded independently, or use the `download_and_rename_genomes.sh` script in the `src\` directory.

1. Open a Julia console and do:
   ```{julia}
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

Next, you can run 

```{bash}
julia --threads <nthreads> scripts/percs.jl -
```