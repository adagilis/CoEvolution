
using Phylo
using PhyloNetworks
using PhyloCoalSimulations
using DataFrames
using Statistics
using HypothesisTests
using CSV

function phylo_to_net(tree)
    return(PhyloNetworks.readnewick(Phylo.outputtree(tree,Newick())))
end

function net_to_phylo(net)
    return(Phylo.parsenewick(PhyloNetworks.writenewick(net)))
end

function run_simulations(species_tree,scale,n)
    copy_tree = deepcopy(species_tree)
    for e in getbranches(copy_tree)
        e.length = scale*e.length
    end
    net = phylo_to_net(copy_tree)
    nets = simulatecoalescent(net,n,1)
    nets = net_to_phylo.(nets)
    for t in nets
        for e in getbranches(t)
            e.length = e.length / scale
        end
    end
    return(nets)
end