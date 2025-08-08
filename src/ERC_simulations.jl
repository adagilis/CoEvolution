
using Phylo
using PhyloNetworks
using PhyloCoalSimulations
using DataFrames
using ProgressMeter
using StatsBase
using Distributions
using LinearAlgebra

include("Phylo_utilities.jl")
include("Stats_utilities.jl")

"""
simulate_tree_set(species_tree,scale,genes,min_tips,max_tips) -> Vector{Trees}

Takes a species tree and simulates a bunch of trees using a random subsample of tree tips for each simulation.
The scale parameter is meant to rescale when the tree is not in coalescent units, but does not seem to work, so keep set to 1.

"""

function simulate_tree_set(species_tree,scale,genes,min_tips,max_tips)
    tips = rand(min_tips:max_tips,genes)
    gene_trees = [simulate_subsampled_tree(species_tree,tips[i],scale) for i in 1:genes]
    return(gene_trees)
end

function simulate_subsampled_tree(species_tree,n,scale)
    copy_tree = deepcopy(species_tree)
    leaves = getleafnames(copy_tree)
    rand_tips = sample(leaves,n,replace=false)
    keeptips!(copy_tree,rand_tips)
    rescale_tree!(copy_tree,scale)
    net = phylo_to_net(copy_tree)
    sim = simulatecoalescent(net,1,1;nodemapping=true)[1]
    tree = net_to_phylo(removedegree2nodes!(sim))
    rescale_tree!(tree,1/scale)
    return(tree)
end

function calculate_ERC_exp_null(species_tree,n;cutoff=5)
    bls = get_branch_lengths(species_tree)
    #random, iid sampling of branches, normalization
    sampled_bls = zeros(length(bls),n)
    for i in 1:length(bls)
        rate = Exponential(bls[i])
        sbls = rand(rate,n)
        sampled_bls[i,:] = sbls./bls[i]
    end
    #calculate ERC for whole dataset
    cor_table = cor_test_matrix(sampled_bls,2;cutoff=cutoff)
    return(cor_table)
end