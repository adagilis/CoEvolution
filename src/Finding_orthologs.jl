using CSV, DataFrames

#general functions to generate diamond calls and identify reciprocal best hits from outputs

function diamond_blastp(db,query,path)
    outname = db*"_"*query*".blastp.tsv"
    dbname = db*".dmnd"
    try
        isfile(dbname) || diamond_makedb(db)
        cmd = `diamond blastp -f 6 --iterate -k 1 -d $dbname -q $query -o $outname`
        run(cmd)
    catch
        """
        Could not run blastp using diamond. Check if database exists, and query files are correctly specified.
        """
    end  
end

function diamond_makedb(db,path)
    loc= path*db
    cmd = `diamond makedb --in $loc -o $path/$db`
    try
        run(cmd)
    catch
        """
        Could not run diamond. Check if sequence exists in working directory!
        """
    end
end

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
        rbh = DataFrame(:species1 = species1,:species2 = species2,:i = complete_cases.Column1,:j = complete_cases.Column2)
        return(rbh)
    catch
        """
        Could not find reciprocal best hits. If no diamond errors were thrown, this may simply mean no reciprocal best hits exist.
        """
    end
end

function build_orth_files(rbh,path)
    uniqseq = unique(rbh.i)
    for i in uniqseq
        subrbh = rbh[rbh.i==i,:]
        infile = path*subrbh.species1[1]
        outfile = path*subrbh.species1[1]*subrbh.i[1]*".fasta"
        #grab sequence from species fasta
        run(`echo ">$i >> $outfile`)
        run(`grep -A1 '$i' $infile >> $outfile`)
        for orth in 1:length(subrbh.i)
            seq_name = subrbh.species2[orth]

            #append sequence to outfile
        end
    end
end