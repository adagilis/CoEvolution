#!/usr/bin/env nextflow
/*
* This module identifies reciprocal best hits for all sets of species. It's really just launching a bunch of diamond comparisons. I think this is nominally faster in our full `Julia` implementation, but there is some convenience for running on clusters found here. 
*/

process RBH_table {

    script:
    """
    #!/usr/bin/env julia
    using DrWatson
    @quickactivate
    include(srcdir("Finding_orthologs.jl"))
    download_and_prep_sequences("${taxon}";level="${level}",accession="${abbrev}",dir="./")
    """
    output:
}

workflow {

    main:

}

output {

}