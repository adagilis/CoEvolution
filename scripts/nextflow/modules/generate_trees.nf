#!/usr/bin/env nextflow
/*
* This module 1) creates sequence alignments for each set of reciprocal best hits, 2) runs iqtree on each tree 3) generates a consensus tree with astral
*/


process runIQTREE{
    input
    path aligned_fasta
    val models

    script:
    """
    iqTree3 -s ${aligned_fasta} -nt AUTO
    """
    output:
       
}



process runAstral{

    script:


    output:
        path astral_consensus.newick
        path astral_consensus_rerooted.newick
}


workflow {

    main:
}

