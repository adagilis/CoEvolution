using CSV, DataFrames
using ProgressMeter
using FLoops
using Term


"""
    download_and_prep_sequences(taxon) -> sequence files in /seqs/ folder
"""
function download_and_prep_sequences(taxon;level="chromosome")
    isdir(data_dir*"seqs/") || mkdir(data_dir*"seqs/")
    run(`datasets download genome taxon $taxon --assembly-level $level --dehydrated --include protein,gff3 --filename $data_dir/$taxon.zip --assembly-source "RefSeq" --mag exclude --exclude-atypical`)
    run(`unzip $data_dir/$taxon.zip -d $data_dir/seqs/`)
    run(`datasets rehydrate --directory $data_dir/seqs/`)
    files = filter(contains("GCF"),readdir(data_dir*"seqs/ncbi_dataset/data/",join=true))
    GCFs = filter(contains("GCF"),readdir(data_dir*"seqs/ncbi_dataset/data/",join=false))
    for g in 1:length(GCFs)
        gcf = GCFs[g]
        faa = files[g]*"/protein.faa"
        gff = files[g]*"/genomic.gff"
        name = split(read(`bash src/get_name.sh $gcf`,String)," ")
        out_faa = data_dir*"seqs/"*name[1][1]*"_"*replace(name[2],r"\n"=>"")*".faa"
        out_gff = data_dir*"seqs/"*name[1][1]*"_"*replace(name[2],r"\n"=>"")*".gff"
        run(`mv $faa $out_faa`)
        run(`mv $gff $out_gff`)
    end
    run(`rm -rf $data_dir/seqs/ncbi_dataset`)
end




#general functions to generate diamond calls and identify reciprocal best hits from outputs

"""
    diamond_blastp(db,query,path) 
creates a command to run a diamond blastp search between a db and query species. Assumed that path contains .faa sequences for the db and query species, generates a diamond database if one does not exist already.
"""
function diamond_blastp(db,query,path)
    outname = path*"blast/"*db*"_"*query*".blastp.tsv"
    dbname = path*"blast/"*db*".dmnd"
    query_seq= path*"seqs/"*query*".faa"
    oldstd=stderr
    try
        redirect_stderr(devnull)
        isfile(dbname) || diamond_makedb(db,path)
        cmd = `diamond blastp --threads 4 -f 6 --iterate -k 1 -d $dbname -q $query_seq -o $outname`
        run(cmd)
        redirect_stderr(oldstd)
    catch
        redirect_stderr(oldstd)
        tprint("{red}Coold not run diamond! Check if properly installed, and sequence exists in directory.{/red}")
    end  
end

"""
    diamond_makedb(db,path) 
creates a command generate a diamond db. Assumed that path contains .faa sequences for the db species
"""
function diamond_makedb(db,path)
    isdir(path*"blast/") || mkdir(path*"blast/")
    loc= path*"seqs/"*db*".faa"
    out= path*"blast/"*db
    cmd = `diamond makedb --in $loc -d $out`
    oldstd=stderr
    try
        redirect_stderr(devnull)
        run(cmd)
        redirect_stderr(oldstd)
    catch
        redirect_stderr(oldstd)
        tprint("{red}Coold not make diamond db! Check if properly installed, sequence exists in directory, and new files can be made.{/red}")
    end
end

"""
    rbh(species1,species2,path) -> DataFrame(:species1,:species2,:gene_species1,:gene_species2)
Identifies reciprocal best hits between species1 and species2, given .faa sequences in `path`. Runs diamond, returns a dataframe
"""
function rbh(species1,species2,path)
    filename1 = path*"blast/"*species1*"_"*species2*".blastp.tsv"
    filename2 = path*"blast/"*species2*"_"*species1*".blastp.tsv"
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
        tprint("{red}Could not detect any reciprocal best hits.{/red}")
    end
end

"""
    find_orths(species1,path) 
Identifies reciprocal best hit orthologs for all genes in species1 versus all sequence files found in path.
Long run times since we allow diamond to take up all threads rather than parallelizing the process here.
Because we are not doing all to all comparisons, the run time is significantly lower than OrthoFinder, however. 
If you want to run this more efficiently, consider running the rbh(species1,species2,path) function across multiple machines/nodes.
The current implementation suffers from having to vcat each result - this is slow, but we don't know a priori how mucn memory to pre-allocate.
"""
function find_orths(focal,species_list,path)
    p=Progress(length(species_list)-1,desc="Running RBH comparisons")
    rbh_res = DataFrame(species1=String[],species2=String[],i=String[],j=String[])
    for s2 in filter(e->e!=focal,species_list)
        rbh_res=vcat(rbh_res,rbh(focal,s2,path))
        next!(p)
    end
    return(rbh_res)
end

"""
    find_orths(rbh_res,path,cutoff=4) 
Accepts a table that lists reciprocal best hit results and creates aligned sequence files for each set of orthologs.
Will not output sequences when there are fewer than `cutoff` (default=4) species with the ortholog, as these will not be useful for generating branch lengths either way.
"""
function build_orth_files(rbh_res,focal,path;cutoff=4,quiet=true)
    outpath=path*focal*"/aligned/"
    #Clean up from previous failed runs, just in case. Deletes all 0 size files
    run(`find $outpath -size 0b -delete`)
    isdir(outpath) || mkdir(outpath)
    focal = rbh_res.species1[1]
    infile = path*"seqs/"*focal*".faa"
    uniqseq = unique(rbh_res.i)
    tprint("{blue}Found $(length(uniqseq)) sets of reciprocal best hits. Generating alignment files.{/blue}")
    p=Progress(length(uniqseq),desc="Aligning files with sufficient sequences (default = 4). Skipping existing files.")
    @floop ThreadedEx() for i in uniqseq
        subrbh = rbh_res[rbh_res.i.==i,:]
        if length(subrbh.species2) >= cutoff
            outfile = outpath*i*".fasta"
            aligned_out = replace(outfile,r".fasta" => s".aligned.fasta")
            if isfile(aligned_out)
                """
                Found existing $aligned_out, skipping sequence! 
                """
            else
                if !quiet 
                    tprint("{red}$i{/red}")
                end
                #grab sequence from species fasta
                append_seq(i,focal,infile,outfile)
                for orth in 1:length(subrbh.i)
                    species_name = subrbh.species2[orth]
                    species_file = path*"seqs/"*species_name*".faa"
                    seq_name = subrbh.j[orth]
                    append_seq(seq_name,species_name,species_file,outfile)
                    #append sequence to outfile
                end
            	run(pipeline(`mafft --anysymbol --quiet --auto $outfile`;stdout=aligned_out))
            	run(`rm $outfile`)
            end
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
