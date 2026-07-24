#!/usr/bin/env nextflow
/*
* This module identifies reciprocal best hits for all sets of species. It's really just launching a bunch of diamond comparisons. I think this is nominally faster in our full `Julia` implementation since rather than run on all of the comparisons we can check which are done and only run remaining. E.g. when using a new species the approach is pretty fast. But let's make it work for nextflow.
*/

params {
    focal: String
    inputDir: Path
}

process makeDBfocal {
    input: 
        path focal
    
    script:
    """
        diamond makedb --in ${focal} -d ${focal}.dmnd
    """
    output:
    tuple path("${focal}"), path("${focal}.dmnd"), val(focal.name)
}

process makeDBspecies {
    input: 
        path species
    
    script:
    """
        diamond makedb --in ${species} -d ${species}.dmnd
    """
    output:
    tuple path("${species}"), path("${species}.dmnd"), val(species.name)
}

process reciprocal_blast {
    input:
        tuple path(focal_seq), path(focal_db), path(species_seq), path(species_db), val(species_pair)
    
    script:
    """
        diamond blastp --threads 4 -f 6 --iterate -k 1 -d $focal_db -q $species_seq -o ${species_pair}_for.tsv
        diamond blastp --threads 4 -f 6 --iterate -k 1 -d $species_db -q $focal_seq -o ${species_pair}_rev.tsv
    """
    output:
    tuple path("${species_pair}_for.tsv"), path("${species_pair}_rev.tsv"), val(species_pair)
}

process reciprocal_best_hits {
    input:
        tuple path(onetotwo), path(twotoone), val(speciespair)

    script:
    /*
    * There's probably a quick bash solution here, and the `Julia` startup time is not worth it, but I'll return to this later.
    */
    """
    #!/usr/bin/env julia
    using DrWatson
    @quickactivate
    using CSV, DataFrames
    one_to_two = CSV.read("${onetotwo}",DataFrame,header=false)
    two_to_one = CSV.read("${twotoone}",DataFrame,header=false)
    rename!(two_to_one,:Column1 => :Column2,:Column2=>:Column1)
    complete_cases = innerjoin(one_to_two,two_to_one,on=[:Column1,:Column2],makeunique=true)
    species1,species2 = split("${speciespair}","-")
    rbh = DataFrame(:species1 => species1,:species2 => species2,:j => complete_cases.Column1,:i => complete_cases.Column2)
    CSV.write("${speciespair}_rbh.csv",rbh)
    """
    output:
    path "${speciespair}_rbh.csv"
}

process concat_files {
    input: 
    path files
    val focal
    script:
    """
    head -1 ${files[0]} > ${focal}_rbh_table.csv
    tail -n +2 ${files} >> ${focal}_rbh_table.csv
    """
    output:
    path "${focal}_rbh_table.csv", emit: rbh
}

process create_seq_files {
    input:
    path rbh_table
    path inputdir
    val focal
    script:
    """
    #!/usr/bin/env julia
    using DrWatson
    @quickactivate
    using CSV, DataFrames
    rbh_res = CSV.read("${rbh_table}",DataFrame)
    uniqseq = unique(rbh_res.i)
    seqfiles = filter(endswith(".faa"),readdir("${inputdir}"))
    focalfile = filter(contains("${focal}"),seqfiles)
    function append_seq(seq,species,species_file,file)
        run(pipeline(`echo ">\$species"`;stdout=file,append=true))
        run(pipeline(pipeline(`seqkit grep -p "\$seq" \$species_file`,`grep -v "\$seq"`);stdout=file,append=true))
    end
    using FLoops, ProgressMeter
    p=Progress(length(uniqseq),desc="Creating sequence files.")
    @floop ThreadedEx() for i in uniqseq
        subrbh = rbh_res[rbh_res.i.==i,:]
        if length(subrbh.species2) >= 4
            outfile = i*".fasta"
            append_seq(i,focalfile,infile,outfile)
            for orth in 1:length(subrbh.i)
                    species_name = subrbh.species2[orth]
                    species_file = filter(contains(species_name),seqfiles)
                    seq_name = subrbh.j[orth]
                    append_seq(seq_name,species_name,species_file,outfile)
            end
        end
        next!(p)
    end
    """
    output:
    path "*.fasta"

}

workflow {

    main:
    focal_file = Channel.fromPath("${params.inputDir}/${params.focal}.faa",checkIfExists:true)
    //Generate a list of all species with data minus focal species.
    targets_ch = Channel.fromPath("${params.inputDir}/*.faa")
                .filter { it.name != "${params.focal}.faa" }

    focal_db_ch = makeDBfocal(focal_file)
    target_db_ch = makeDBspecies(targets_ch)
    pairs_ch = focal_db_ch.combine(target_db_ch)
            .map { focal_seq, focal_db, focal_id, species_seq, species_db, species_id -> 
                def focal = focal_id.replace(".faa","")
                def target = species_id.replace(".faa","")
                def pair_id = "${focal}-${target}"
                return tuple(focal_seq, focal_db, species_seq, species_db, pair_id)
            }
    
    rb =reciprocal_blast(pairs_ch)
    rbh_ch = reciprocal_best_hits(rb)
    rbh_res = concat_files(rbh_ch.collect(),params.focal)

    //next we use `Julia` to make unaligned fastas of each resulting group of orthologs

    prep_orth_files = create_seq_files(rbh_res.rbh,params.inputDir,params.focal)

    //aligned_files = align_fastas(prep_orth_files.)

    publish:
    rbh = rbh_res.rbh

}

output {
    rbh {
        path 'data/'
    }
}