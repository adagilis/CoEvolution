using CSV, DataFrames
using ProgressMeter

#general functions to generate diamond calls and identify reciprocal best hits from outputs

"""
    diamond_blastp(db,query,path) 
    creates a command to run a diamond blastp search between a db and query species. Assumed that path contains .faa sequences for the db and query species, generates a diamond database if one does not exist already.
"""

function diamond_blastp(db,query,path)
    outname = path*db*"_"*query*".blastp.tsv"
    dbname = path*db*".dmnd"
    query_seq= path*query*".faa"
    try
        isfile(dbname) || diamond_makedb(db,path)
        cmd = `diamond blastp -f 6 --iterate -k 1 -d $dbname -q $query_seq -o $outname`
        run(cmd)
    catch
        """
        Could not run blastp using diamond. Check if database exists, and query files are correctly specified.
        """
    end  
end

"""
    diamond_makedb(db,path) 
    creates a command generate a diamond db. Assumed that path contains .faa sequences for the db species
"""


function diamond_makedb(db,path)
    loc= path*db*".faa"
    out= path*db
    cmd = `diamond makedb --in $loc -d $out`
    try
        run(cmd)
    catch
        """
        Could not run diamond. Check if sequence exists in working directory!
        """
    end
end

"""
    rbh(species1,species2,path) -> DataFrame(:species1,:species2,:gene_species1,:gene_species2)
    Identifies reciprocal best hits between species1 and species2, given .faa sequences in `path`. Runs diamond, returns a dataframe
"""

function rbh(species1,species2,path)
    filename1 = path*species1*"_"*species2*".blastp.tsv"
    filename2 = path*species2*"_"*species1*".blastp.tsv"
    try
        isfile(filename1) || diamond_blastp(species1,species2,path)
        isfile(filename2) || diamond_blastp(species2,species1,path)
        one_to_two = CSV.read(filename1,DataFrame,header=false)
        two_to_one = CSV.read(filename2,DataFrame,header=false)
        rename!(two_to_one,:Column1 => :Column2,:Column2=>:Column1)
        complete_cases = dropmissing(leftjoin(one_to_two,two_to_one,on=[:Column1,:Column2],makeunique=true))
        rbh = DataFrame(:species1 => species1,:species2 => species2,:j => complete_cases.Column1,:i => complete_cases.Column2)
        return(rbh)
    catch
        """
        Could not find reciprocal best hits. If no diamond errors were thrown, this may simply mean no reciprocal best hits exist.
        """
    end
end

"""
    find_orths(species1,path) 
    Identifies reciprocal best hit orthologs for all genes in species1 versus all sequence files found in path.
    Long run times since we allow diamond to take up all threads rather than parallelizing the process here.
    Because we are not doing all to all comparisons, the run time is significantly lower than OrthoFinder, however. 
    If you want to run this more efficiently, consider running the rbh(species1,species2,path) function across multiple machines/nodes.
"""

function find_orths(focal,species_list,path)
    p=Progress(length(species_list)-1,desc="Running RBH comparisons")
    rbh_res = DataFrame(species1=String[],species2=String[],i=String[],j=String[])
    oldstd=stderr
    redirect_stderr(devnull)
    for s2 in filter(e->e!=focal,species_list)
        rbh_res=vcat(rbh_res,rbh(focal,s2,path))
        next!(p)
    end
    redirect_stderr(oldstd)
    return(rbh_res)
end

"""
    find_orths(rbh_res,path,cutoff=4) 
    Accepts a table that lists reciprocal best hit results and creates aligned sequence files for each set of orthologs.
    Will not output sequences when there are fewer than `cutoff` (default=4) species with the ortholog, as these will not be useful for generating branch lengths either way.
"""

function build_orth_files(rbh_res,path;cutoff=3)
    outpath=path*"aligned/"
    isdir(outpath) || mkdir(outpath)
    focal = rbh_res.species1[1]
    infile = path*focal*".faa"
    uniqseq = unique(rbh_res.i)
    p=Progress(length(uniqseq),desc="Aligning files with sufficient sequences")
    for i in uniqseq
        subrbh = rbh_res[rbh_res.i.==i,:]
        if length(subrbh.species2) >= cutoff
            outfile = outpath*i*".fasta"
            #grab sequence from species fasta
            append_seq(i,focal,infile,outfile)
            for orth in 1:length(subrbh.i)
                species_name = subrbh.species2[orth]
                species_file = path*species_name*".faa"
                seq_name = subrbh.j[orth]
                append_seq(seq_name,species_name,species_file,outfile)
                #append sequence to outfile
            end
            aligned_out = replace(outfile,r".fasta" => s".aligned.fasta")
            run(`muscle -align $outfile -output $aligned_out`)
            run(`rm $outfile`)
        end
        next!(p)
    end
end

"""
    append_seq(seq,species,species_file,file) 
    Appends the sequence `seq` from `species` into `file`. Requires `species_file` to exist, and needs `seqkit`.
    Appends the sequence as: 
    ```{fasta}
    >species
    Sequence
    ```
"""

function append_seq(seq,species,species_file,file)
    run(pipeline(`echo ">$species"`;stdout=file,append=true))
    run(pipeline(pipeline(`seqkit grep -p "$seq" $species_file`,`grep -v "$seq"`);stdout=file,append=true))
end