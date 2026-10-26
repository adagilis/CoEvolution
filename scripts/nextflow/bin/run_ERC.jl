using DrWatson
@quickactivate
using ArgParse

function parse_cmd()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--inputDir", "-d"
            help= "Directory containing tree files"
            required=true
        "--focal","-f"
            help="Focal species id"
            default=datadir()
            required=true
    end
    return parse_args(s)
end

parsed_args = parse_cmd()
focal = parsed_args["focal"]
dir = parsed_args["inputDir"]
include(srcdir("ERC_functions.jl"))

    species_tree = read_tree(dir*"/"*focal*"_rerooted.newick.clean")
    tree_files = filter(endswith(".treefile"),readdir(dir,join=true))

    #This takes a while, but saves a lot of reading time in the ERC loop. It does increase memory requirements even further, consider running `runERC_files(trefiles,species_tree)` if you are being memory limited. But, you are avoiding sum(1:length(tree_files)) file read operations, which is also probably a benefit.
    trees = [read_tree(t) for t in tree_files]
    ERC = runERC_collection(trees,species_tree;cutoff=${params.cutoff},min_shared=${params.minShared})
    ERC = ERC[findall(completecases(ERC)),:]
    disallowmissing!(ERC)
    genes = replace.(filter(endswith(".treefile"),readdir(dir,join=false)),r"aligned.fasta.treefile"=>s"")

    ERC.i = genes[Int.(ERC.i)]
    ERC.j = genes[Int.(ERC.j)]
    ERC.gid = [join(sort([x.i,x.j]),"_") for x in eachrow(ERC)]

    using Arrow
    Arrow.write(focal*"_ERC.arrow",ERC)