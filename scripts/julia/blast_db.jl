using DrWatson
@quickactivate
using ArgParse

function parse_cmd()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--focal", "-f"
            help= "Either an NCBI taxon id (recommended), or a taxon name. Note that even some seemingly obvious taxon names may have duplicates in other groups, so use taxon IDs when possible."
            required=true
        "--output_dir","-o"
            help="Output directory. By default, data is placed in: './data/focal/'." 
            default=""
            required=false
        "--seqs_dir","-s"
            help="Directory containing .faa files of species to scan." 
            default=datadir("seqs/")
            required=true
    end
    return parse_args(s)
end

parsed_args = parse_cmd()

"""
check_for_diamond() -> Boolean
Quick function to check if datasets is installed and functioning.
"""
function check_for_diamond()
    test = run(`diamond --version`)
    if(test.exitcode==0)
        return(true)
    else
        return(false)
    end
end


check_for_diamond() || """`diamond` could not be run, check if it is installed"""

data_dir = parsed_args["dir"]
if data_dir==""
    data_dir=datadir(parsed_args["focal"]*"/")
elseif( !endswith(data_dir,"/"))
    data_dir *= "/"
end
species = replace.(filter(endswith(".faa"),readdir(parsed_args["seqs_dir"]))),r".faa"=>"")
all_orthologs = find_orths(parsed_args["focal"],species,parsed_args["seqs_dir"],data_dir)

using Arrow

Arrow.write(data_dir*focal*"_orthologs.arrow",all_orthologs)

