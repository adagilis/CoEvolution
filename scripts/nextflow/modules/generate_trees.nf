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
    iqtree3 -s ${aligned_fasta} -nt AUTO
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
    astral4 -i concat_trees.treefile -o ${params.focal}_consensus.newick --root ${outgroup}
    """
    output:
        path "${params.focal}_consensus.newick", emit: tree
}

process reRootTree{
    input:
    path tree
    val focal

    script:
    """
    #!/usr/bin/env Rscript
    library(ape)
    tree = read.tree("${tree}")
    rerooted = root.phylo(tree,outgroup="${focal}")
    write.tree(rerooted,file="${focal}_rerooted.newick")
    """
    output:
        path "${focal}_rerooted.newick", emit: tree
}

process cleanTree{
    input:
    path tree

    script:
    """
    sed -i 's/NaN/0.0/g' ${tree}
    sed  's/)[0-1].[0-9]*:/):/g' ${tree} > ${tree}.clean
    """
    output:
    path "*_rerooted.newick.clean",emit:tree
}

workflow {

    main:
    aligned_genes = Channel.fromPath("${params.inputDir}/*.aligned.fasta")
    make_phylo = runIQTREE(aligned_genes)

    make_consensus = runAstral(make_phylo.tree.collect(),params.outgroup)
    rerooted = reRootTree(make_consensus.tree,params.focal)
    cleaned = cleanTree(rerooted.tree)

    publish: 
    trees = make_phylo.tree
    consensus = make_consensus.tree
    rerooted = cleaned.tree
}

output {
    trees{
        path "data/trees/"
    }
    consensus{
        path "data/trees/"
    }
    rerooted{
        path "data/trees/"
    }
}