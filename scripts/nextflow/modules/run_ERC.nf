#!/usr/bin/env nextflow
/*
* This module takes a set of gene trees and a species tree and returns ERC scores in an Arrow file
*/

params {
    focal: String
    inputDir: Path
    cutoff: Integer = 5
    minShared: Integer = 4
}

process runERC{
    cpus 12
    input:
        path inputDir
        val focal
    script:
    """
    julia --threads 12 run_ERC.jl --focal "${params.focal}" --inputDir "${params.inputDir}"
    
    """

    output:
    path "*.Arrow"
}

workflow {

    main:
    ERC = runERC(params.inputDir,params.focal)


    publish:
    ERC_res = ERC.out
}

output{
    ERC_res{
        path 'data/'
    }
}
