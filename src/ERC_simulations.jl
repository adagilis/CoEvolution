
using Phylo
using PhyloNetworks
using PhyloCoalSimulations
using DataFrames
using ProgressMeter
using StatsBase

include("Phylo_utilities.jl")

function run_simulations(species_tree,scale,genes,min_tips,max_tips)
    tips = rand(min_tips:max_tips,genes)
    gene_trees = [simulate_subsampled(species_tree,tips[i],scale) for i in 1:genes]
    return(gene_trees)
end

function simulate_subsampled(species_tree,n,scale)
    local copy_tree = deepcopy(species_tree)
    local leaves = getleafnames(copy_tree)
    local rand_tips = sample(leaves,n,replace=false)
    reorder_tree(copy_tree,rand_tips)
    rescale_tree!(copy_tree,scale)
    local net = phylo_to_net(copy_tree)
    local sim = simulatecoalescent(net,1,1)[1]
    local tree = net_to_phylo(sim)
    rescale_tree!(tree,1/scale)
    return(tree)
end