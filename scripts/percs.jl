using DrWatson
@quickactivate "CoEvolution"

"""
Percs v0.0.1
#I hate this name, but I need version numbers and feels like a version number to an untitled piece of software is horrible.

Pipeline for Evolutionary Rate Correlation Searches

This software will calculate ERC values for genes from a target species to all species with supplied protein sequences.

    Requires DIAMOND, iqTree2

Step 1: Initializing packages and precompiling code.

"""

using ArgParse, JLD2
function parse_arguments()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--focal_species","-f"
            help = "Focal species name. (should exist in protein sequence folder supplied)"
            required=true
            arg_type = String
        "--outgroup_species","-o"
            help = "Outgroup species name. (should exist in protein sequence folder supplied)"
            required=true
            arg_type = String
        "--species_tree", "-t"
            required=true
            help = "Newick file containing species tree of relevant species if new one is not desired."
        "--seq_dir", "-s" 
            required = true
            help = "Directory containing protein sequences for each species of interest. Files should have suffix of *.faa or *.fasta"
            arg_type = String
        "--redo"
            help = "Flag to ignore existing outputs and generate new ones."
            action = :store_true
        "--rescale", "-r"
            help = "Flag if branch lengths need to be rescaled, especially if supplying own species tree."
            action = :store_true
        "--theta", "-θ"
            arg_type = Float64
            default = 1.0
            help = "Watterson's θ used to rescale branch lengths when needed (default = 1)."
        "--simulations", "-s"
            help = "number of gene trees to generate under exponential assumptions (default = 1000)"
            arg_type = Int
            default = 1000
        "--project_name","-p"
            required = true
            help = "Folder to output ERC results, simulations and plots."
    end
    return parse_args(s)
end

parsed_args = parse_arguments()

isdir(parsed_args["project_name"]) || mkdir(parsed_args["project_name"])

include(srcdir("ERC_functions.jl"))
include(srcdir("ERC_simulations.jl"))
include(srcdir("Graphs_utilities.jl"))

"""
Successfuly precompiled!

Currently active project is: $(parsed_args["project_name"])

Path of active project: $(projectdir())

Outputs being written to: $(parsed_args["project_name"])

Examining project folder.
"""

#Look for existing files, and skip to step depending on what exists unless there's a redo flag.

#checkpoint states
states = Dict(
    :current_step => 1,
    :protein_fastas_found => Vector{String},
    :outgroup_files_found => Vector{String},
    :focal_files_found => "",
    :orthologs_detected => 0,
    :trees_found => 0,
)

#Look for checkpoint
ckp_file = parsed_args["project_name"]*"/runtime.ckp"
isfile(ckp_file) && states = jldopen(ckp_file)["states"]

while states[:current_step] == 1
    """
    Could not find previous run. Searching for all relevant files.
    """
    species = filter(endswith(".faa") || contains(".fasta"),readdir(parsed_args["seq_dir"],join=true))
    states[:protein_fastas_found] = species
    if stats[:protein_fastas_found] == 0
        """
        Could not detect protein fastas. Check the inputs!
        """
        break
    end

    focal = filter(contains(parsed_args["focal_species"]),species)
    outgroup = filter(contains(parsed_args["outgroup_species"]),species)

    if  length(focal) == 0
        """
        Could not find focal species. Check for typos?
        """
        break
    elseif length(focal) > 1
        """
        Found multiple protein files matching supplied focal species name. Cannot proceed due to ambiguous focal species.
        """
        break
    elseif length(focal) == 1
        """
        Found focal species.
        """
        states[:focal_files_found] = true
    end

    if  length(outgroup) == 0
        """
        Could not find outgroup species. Check for typos?
        """
        break
    elseif length(outgroup) > 1
        """
        Found multiple protein files matching supplied outgroup species name. This has not been formally tested!
        """
        states[:outgroup_files_found] = outgroup
        states[:current_step] = 2
    elseif length(outgroup) == 1
        """
        Found outgroup species.
        """
        states[:outgroup_files_found] = outgroup
        states[:current_step] = 2
    end
end

#Step 2: Run Diamond.

while states[:current_step]==2
    include(srcdir("finding_orthologs.jl"))


    states[:current_step] = 3
end

#Create trees, starting with species tree.
while states[:current_step] == 3
    #Requires iqTree2 to be callable.
    try

    catch
    end
end

#Calculate ERC scores
while states[:current_step] == 4
    """
    Running step 4: calculating ERC scores
    """
    
    trees = filter(contains(".treefile"),readdir(parsed_args["tree_dir"],join=true))
    species_tree = read_tree(parsed_args["species_tree"])

    println(
    """
    Found $(length(trees)) in $(parsed_args["tree_dir"]), loaded in species tree.

    Running $(binomial(length(trees),2)) comparisons.
    """
    )


    ERC = runERC_files(trees,species_tree)

    #Quick report of results
    using UnicodePlots

    println("""
    ERC values calculated! μ = $(mean(ERC[:,"r"])), σ = $(std(ERC[:,"r"])).

    $(length(findall(ERC.n_edges .< 4))) interactions excluded due to too few edges in trees overlapping.

    Distribution:
    """)

    histogram(ERC[ERC.n_edges .> 0,"r"])

    #Save output into jld2

    using JLD2

    jldsave(parsed_args["project_name"]*"/ERC_stats.jld2",ERC=ERC)
    states[:current_step] = 5
end

while states[:current_step] ==5
    """
    Project saved! 

    Step 5: running null simulations based on species_tree for $(parsed_args["simulations"]) simulations
    """

    #Determine cutoffs based on null

    ERC_null = calculate_ERC_exp_null(species_tree,parsed_args["simulations"])

    jldsave(parsed_args["project_name"]*"/ERC_nulls.jld2",ERC_null=ERC_null)

    limits = quantile(ERC_null.r,[0.0275,0.975])

    """
    Null simulations produced, cut-offs for top/bottom 2.5% are: $(limits).

    Discarding values withing 95% of null in real data results in $(length(intersect(findall(ERC.r .> limits[1]),findall(ERC.r .< limits[2])))) values.
    """

    sub_ERC = ERC[union(findall(ERC.r .< limits[1]),findall(ERC.r .> limits[2])),:]

    states[:current_step] = 6
end
#Generating a network

"""
Step 4: Generating network. 

Assuming genes interact if outside 95% expected from null, and p-val < 0.05 / $(length(sub_ERC.r)).

Results in $(length(findall(sub_ERC.pval .< 0.05/length(sub_ERC.r)))) edges.

Significant edges output to: $(parsed_args["project_name"]*"/network.tsv"))
"""

sub_ERC = ERC[ERC.pval .< 0.05/length(sub_ERC.r),:]

CSV.write(parsed_args["project_name"]*"/significant_edges.csv",sub_ERC)

#include(srcdir("Graphs_utilities.jl"))

#graph = construct_ERC_graph(sub_ERC)

#writecsv

