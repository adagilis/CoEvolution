
using Phylo
using PhyloNetworks
using PhyloCoalSimulations
using DataFrames
using Statistics
using HypothesisTests
using CSV

include("Phylo_utilities.jl")

function run_simulations(species_tree,scale,n)
    copy_tree = deepcopy(species_tree)
    rescale_tree!(copy_tree,scale)
    net = phylo_to_net(copy_tree)
    nets = simulatecoalescent(net,n,1)
    nets = net_to_phylo.(nets)
    for t in nets
        rescale_tree!(t,1/scale)
    end
    return(nets)
end