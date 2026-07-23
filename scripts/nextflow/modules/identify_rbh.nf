#!/usr/bin/env nextflow
/*
* This module identifies reciprocal best hits for all sets of species. It's really just launching a bunch of diamond comparisons. I think this is nominally faster in our full `Julia` implementation since rather than run on all of the comparisons we can check which are done and only run remaining. E.g. when using a new species the approach is pretty fast. But let's make it work for nextflow.
*/

params {
    focal: String
    inputDir: Path
}

process makeDB {
    input: 
        path species
    
    script:
    def out_name = species.replace(".faa",".dmnd")
    """
        diamond makedb --in ${species} -d ${out_name}
    """

    output:
    path "*.dmnd"
}

process reciprocal_blast {
    input:
        path species1_seq
        path species2_seq
        path species1_db
        path species2_db
    
    script:
    """
        diamond blastp --threads 4 -f 6 --iterate -k 1 -d $species1_db -q $species2_seq -o 1_2.tsv
        diamond blastp --threads 4 -f 6 --iterate -k 1 -d $species2_db -q $species1_seq -o 2_1.tsv
    """
    output:
    path "2_1.tsv", emit: 2_1
    path "1_2.tsv", emit: 1_2
}

process reciprocal_best_hits {
    input:
        path 1_2
        path 2_1
    
    script:
    /*
    * There's probably a quick bash solution here, and the `Julia` startup time is not worth it, but I'll return to this later.
    */
    """
    #!/usr/bin/env julia
    using DrWatson
    @quickactivate
    using CSV, DataFrames
    one_to_two = CSV.read(filename1,DataFrame,header=false)
    two_to_one = CSV.read(filename2,DataFrame,header=false)
    rename!(two_to_one,:Column1 => :Column2,:Column2=>:Column1)
    complete_cases = innerjoin(one_to_two,two_to_one,on=[:Column1,:Column2],makeunique=true)
    rbh = DataFrame(:species1 => species1,:species2 => species2,:j => complete_cases.Column1,:i => complete_cases.Column2)
    CSV.write("rbh.csv",rbh)
    """
    output:
    rbh.csv
}

process RBH_table {
    input: 
    path 'csv'

    script:
    """
    #!/usr/bin/env julia
    using DrWatson
    @quickactivate
    using CSV, DataFrames, Arrow
    csvs = filter(endswith(".csv"),readdir("."))
    rbh_res = DataFrame(species1=String[],species2=String[],i=String[],j=String[])
    for x in 1:length(csvs)
        tmp = CSV.read(csvs[x],DataFrame)
        rbh_res=vcat(rbh_res,tmp)
    end
    Arrow.write("orthologs.Arrow",rbh_res)
    """
    output:
    path "orthologs.Arrow"
}

workflow {

    main:
    focal_file = file("${params.inputDir}/${params.focal}.faa")
    //Generate a list of all species with data minus focal species.
    targets_ch = Channel.fromPath("${params.inputDir}/*.faa")
        .filter { it.name != focal_file.name }

    focal_db_ch = makeDB(focal_file)
    target_db_ch = makeDB(targets_ch)
    rbh_ch = focal_file.combine(targets_ch)
            .combine(focal_db_ch)
            .combine(target_db_ch)




}

output {

}