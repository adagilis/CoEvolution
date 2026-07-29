#!/usr/bin/env nextflow
/*
* This module 1) creates sequence alignments for each set of reciprocal best hits, 2) runs iqtree on each tree 3) generates a consensus tree with astral
*/

params {
    focal: String
    inputDir: Path
    outgroup: String
}

process runIQTREE{
    cpus 4

    input:
    path aligned_fasta

    script:
    """
    iqTree3 -s ${aligned_fasta} -nt AUTO
    """
    output:
    path "*.treefile",emit:tree
}



process runAstral{
    cpus 8

    input:
    path trees
    val outgroup

    script:
    """
    cat ${trees} > concat_trees.treefile
    astral4 -i concat_trees.treefile -o astral_consensus.newick --root ${outgroup}

    """
    output:
        path astral_consensus.newick, emit: consensus
        path astral_consensus_rerooted.newick, emit: rerooted
}


workflow {

    main:
    aligned_genes = Channel.fromPath("${params.inputDir}/*.aligned.fasta")
    make_phylo = runIQTREE(aligned_genes)

    make_consensus = runAstral(make_phylo.tree.collect())

    publish: 
    make_phylo.tree
    make_consensus.consensus
    make_consesnsus.rerooted
}

