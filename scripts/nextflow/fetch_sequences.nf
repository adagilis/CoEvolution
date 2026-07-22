#!/usr/bin/env nextflow
/*
*
*/

process fetchSequences {
    script:
    """
    julia scripts/julia/fetchSequences.jl -s $species -l $level -d $path -a $name
    """

    output:
        path 'seqs/*.faa'
}

workflow {

    main:
    // emit a greeting
    fetchSequences()

    publish:
    first_output =  
}