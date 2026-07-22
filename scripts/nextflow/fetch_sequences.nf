#!/usr/bin/env nextflow
/*
*
*/

params {
    level: String = "chromosome"
    abbrev: String = "classic"
}

process fetchSequences {
    input:
    val taxon
    val level
    val abbrev

    script:
    """
    julia scripts/julia/fetchSequences.jl -t '${taxon}' -l '${level}' -d "." -a '${abbrev}'
    """

    output:
    path "seqs/*.faa"
}

workflow {

    main:
    // emit a greeting
    fetchSequences(params.taxon,params.level,params.abbrev)

}