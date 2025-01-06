using DrWatson
@quickactivate "CoEvolution"

# Here you may include files from the source directory
include(srcdir("ERC_functions.jl"))

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())

Running analysis step 1: calculating ERC scores
"""
)

trees = filter(contains(".treefile"),readdir(datadir("trees"),join=true))
species_tree = open(parsenewick,datadir("trees","species_tree.newick"))

println(
"""
Found $(length(trees)) in `data/trees` directory, loaded in species tree.

Running $(length(trees)*(length(trees)-1)/2) comparisons.
"""
)

ERC_res = runERC()

