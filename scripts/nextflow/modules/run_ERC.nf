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
    input:
        path inputDir
        val focal
    
    script:
    """
    #!/usr/bin/env julia
    using DrWatson
    @quickactivate
    include(srcdir("ERC_functions.jl"))

    species_tree = read_tree("${inputDir}/${focal}_rerooted.newick.clean")
    tree_files = filter(endswith(".treefile"),readdir("${inputDir}/",join=true))

    #This takes a while, but saves a lot of reading time in the ERC loop. It does increase memory requirements even further, consider running `runERC_files(trefiles,species_tree)` if you are being memory limited. But, you are avoiding sum(1:length(tree_files)) file read operations, which is also probably a benefit.
    trees = [read_tree(t) for t in tree_files]
    ERC = runERC_collection(trees,species_tree;cutoff=${params.cutoff},min_shared=${params.minShared})
    ERC = ERC[findall(completecases(ERC)),:]
    disallowmissing!(ERC)
    genes = replace.(filter(endswith(".treefile"),readdir("${inputDir}/",join=false)),r"aligned.fasta.treefile"=>s"")

    ERC.i = genes[Int.(ERC.i)]
    ERC.j = genes[Int.(ERC.j)]
    ERC.gid = [join(sort([x.i,x.j]),"_") for x in eachrow(ERC)]

    using Arrow
    Arrow.write("${focal}_ERC.arrow",ERC)
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
