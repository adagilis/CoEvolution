using DrWatson
@quickactivate "CoEvolution"

using ArgParse, JLD2
function parse_arguments()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--species_tree", "-t"
            required = true
            help = "Newick file containing species tree to simulate along."
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
        "--positive_sim_branches","-n"
            help = "Number of branches with correlated rates for positive simulation subset (default 4)."
            default=4
            arg_type = Int
    end

    return parse_args(s)
end


#= This script takes a species tree (with branch lengths in coalescent units) 
and simulates a set of ERC metrics for a bunch of genes evolving neutrally 
along the species tree. Then, a set of genes with correlated changes to evolutionary rates along _n_ branches
are simulated. 

The resulting distributions of ERC scores can be used as a  baseline for ERC scores calculated for the actual
data.
=#

include(srcdir("ERC_functions.jl"))
include(srcdir("ERC_simulations.jl"))
include(srcdir("Phylo_utilities.jl"))

parsed_args = parse_arguments()
println("Parsed args:")
for(arg,val) in parsed_args
    println(" $arg => $val")
end

println(
"""
Null simulations for ERC.

Currently active project is: $(projectname())

Reading in tree at $(parsed_args["species_tree"])
""")

tree = read_tree(parsed_args["species_tree"])

#Rescale branch lengths if needed

if(parsed_args["rescale"])
    println(
        """
        Running analysis step: Convert branch lengths to coalescent units.

        Assuming branch lenghts are substitutions/site, and that

        Coalescent units = substitution rate/θ.

        θ should be provided as a parameter, can be calculated as 4 Ne μ.
        """
    )
    rescale_tree!(tree,1/parsed_args["theta"])
end

#Run simulations

println(
"""
Running analysis step: Simulate neutral gene trees.

Simulating $(parsed_args["simulations"]) trees.
""")

trees_simulated = run_simulations(tree,1,parsed_args["simulations"])
println("Trees simulated, calculating $(binomial(parsed_args["simulations"],2)) ERC values.")
ERC_simulated = runERC_collection(trees_simulated,tree)

using UnicodePlots
println("ERC calculated, value distribution:")

println("$(histogram(ERC_simulated[:,"r2"]))")

println("Saving simulations to $(datadir("sims"))")

jldsave(datadir("sims","simulation_results.jld2"); trees=trees_simulated,ERC=ERC_simulated)

println("")

#Simulate small number of trees with correlations at n branches.



println("Done!")