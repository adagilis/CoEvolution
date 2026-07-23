#!/usr/bin/env nextflow
/*
* This module takes an ncbi taxon ID and downloads all amino acid sequences of chromosome level assemblies for that taxon. It then runs `primary_transcript.py` to reduce the data to primary transcripts.def () {
* We run this module to download both the focal species at any level of assembly, outgroup and high quality species within the taxon more broadly. No current module to identify a good set of species more broadly   
}
*/
params {
    taxon: String
    outdir: Path = "."
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
    download_and_prep_sequences("${taxon}";level="${level}",accession="${abbrev}",dir="./")
    """

    output:
    path "*.faa", emit: faas
    path "*.gff", emit: gffs
}

process primaryTranscripts {

    input:
    path species

    script:
    """
    python3 ~/OrthoFinder_source/tools/primary_transcript.py ${species}
    mv primary_transcripts/*.faa ${species}
    rm -rf primary_transcripts
    """

    output:
    path "${species}"
}


workflow {

    main:

    // Grab sequences for taxon and outgroup, so make a channel
    species_ch = fetchSequences(params.taxon,params.level,params.abbrev)

    //Find primary transcripts for each resulting .faa
    primary_ch = primaryTranscripts(species_ch.faas.flatten())

    publish:
    primary = primary_ch
    gffs = species_ch.gffs

}

output {
    primary {
        path 'data/seqs'
    }
    gffs {
        path 'data/seqs'
    }
}