using DrWatson
@quickactivate
using ArgParse

function parse_cmd()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--taxon", "-t"
            help= "Either an NCBI taxon id (recommended), or a taxon name. Note that even some seemingly obvious taxon names may have duplicates in other groups, so use taxon IDs when possible."
            required=true
        "--dir","-d"
            help="Output directory. By default, data is placed in: './data/seqs/'." #Rewrite to be more generic, but now we have clashes between DrWatson and nextflow. Could stop using DrWatson!
            default=datadir()
            required=false
        "--level", "-l"
            help="Optional argument to limit what level of assembly is included. By default, chromosome level assemblies are queried to save on computational time while maximizing quality."
            default="chromosome"
            required=false
        "--accession","-a"
            help="""How should accession names be treated? By default (option 'classic'), accessions are converted to human readable species names. 
            E.g. GCF_000001405.40 is turned into H_sapiens.
            Options are:\n
            accession: use accession as is (GCF_000001405.40);\n
            classic [default](H_sapiens);\n
            full (Homo_sapiens)"""
            default="classic"
            required=false
    end
    return parse_args(s)
end

parsed_args = parse_cmd()

"""
check_for_ncbi_datasets() -> Boolean
Quick function to check if datasets is installed and functioning.
"""
function check_for_ncbi_datasets()
    test = run(`datasets --version`)
    if(test.exitcode==0)
        return(true)
    else
        return(false)
    end
end

check_for_ncbi_datasets() || """NCBI `datasets` could not be run, check if it is installed"""
include(srcdir("Finding_orthologs.jl"))
data_dir = parsed_args["dir"]
if( !endswith(data_dir,"/"))
    data_dir *= "/"
end
download_and_prep_sequences(parsed_args["taxon"];level=parsed_args["level"],accession=parsed_args["accession"],dir=data_dir)
