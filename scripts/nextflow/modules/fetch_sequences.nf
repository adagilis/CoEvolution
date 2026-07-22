#!/usr/bin/env nextflow
/*
*
*/

params {
    taxon: String
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
    #!/usr/bin/env julia
    using DrWatson
    @quickactivate
    include(srcdir("Finding_orthologs.jl"))
    dir=run(`\${PWD}`)
    if( !endswith(dir,"/"))
        dir *= "/"
    end
    download_and_prep_sequences("${taxon}";level="${level}",accession="${abbrev}",dir=dir)
    """

    output:
    path "seqs/*.faa"
}

workflow {

    main:
    // emit a greeting
    fetchSequences(params.taxon,params.level,params.abbrev)

}