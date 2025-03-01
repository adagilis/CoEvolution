using DrWatson
@quickactivate "CoEvolution"

println("Analysis pipeline for Evolutionary Rate Correlations including null simulations.")

using ArgParse, JLD2
function parse_arguments()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--species_tree", "-t"
            required = true
            help = "Newick file containing species tree to sim."
        "--tree_dir", "-d" #TODO rewrite as a concatenated treefile?
            required = true
            help = "Directory containing gene tree files."
        "--rescale", "-r"
            help = "Flag if branch lengths need to be rescaled"
            action = :store_true
        "--theta", "-θ"
            arg_type = Float64
            default = 1.0
            help = "Watterson's θ used to rescale branch lengths when needed (default = 1)"
        "--simulations", "-s"
            help = "number of gene trees to simulate, exponential run time with increased simulations (default = 1000)"
            arg_type = Int
            default = 1000
        "--project_name","-o"
            help = "Folder to output ERC results, simulations and plots."
    end

    return parse_args(s)
end

parsed_args = parse_arguments()

include(srcdir("ERC_functions.jl"))

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())

Running analysis step 1: calculating ERC scores
"""
)

trees = filter(contains(".treefile"),readdir(parsed_args["tree_dir"]),join=true)
species_tree = read_tree(parsed_args["species_tree"])

println(
"""
Found $(length(trees)) in `data/trees` directory, loaded in species tree.

Running $(binomial(length(trees),2)) comparisons.
"""
)


ERC_res = runERC_files(trees,species_tree)


#Quick report of results
using UnicodePlots

println("""
ERC values calculated! μ = $(mean(ERC_res[:,"r2"])), σ = $(std(ERC_res[:,"r2"])).

$(length(findall(ERC_res.n_edges .< 4))) interactions excluded due to too few edges in trees overlapping.

Distribution:
""")

histogram(ERC_res[ERC_res.n_edges .> 0,"r2"])

#Save output into jld2

using JLD2

jldsave(datadir(parsed_args["project_name"],"ERC_stats.jld2"),ERC=ERC_res)

"""
Project saved, now running null simulations based on species_tree for $(parsed_args["simulations"]) simulations
"""

include(srcdir("ERC_simulations.jl"))
#Determine cutoffs based on null

ERC_null = calculate_ERC_exp_null(species_tree,parsed_args["simulations"])

jldsave(datadir(parsed_args["project_name"],"ERC_stats.jld2"),ERC_null=ERC_null)

limits = quantile(ERC_null.r2,[0.0275,0.975])

"""
Null simulations produced, cut-offs for top/bottom 2.5% are: $(limits).

Discarding values withing 95% of null in real data results in $(length(intersect(findall(ERC.r2 .> limits[1]),findall(ERC.r2 .< limits[2])))) values.
"""

sub_ERC = ERC[intersect(findall(ERC.r2 .> limits[1]),findall(ERC.r2 .< limits[2])),:]

#Generating a network

"""
Generating network. Assuming genes interact if p-val < 0.05 / $(length(sub_ERC.r2)).

Results in $(length(findall(ERC_res[:,"pval"] .< 0.05/length(sub_ERC.r2)))) edges.

Significant edges output to: $(datadir("processed","network.tsv"))
"""

sub_ERC = ERC_res[ERC_res.pval .< 0.05/length(sub_ERC.r2),:]

include(srcdir("Graphs_utilities.jl"))

graph = construct_ERC_graph(sub_ERC)

#writecsv





