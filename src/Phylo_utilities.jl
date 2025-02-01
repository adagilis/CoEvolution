using Phylo
using PhyloNetworks

"""
    read_tree(file) → Phylo
Quick utility to read a newick file in as a Phylo object.
"""
function read_tree(file)
    return(open(parsenewick,file))
end

"""
    phylo_to_net(tree:Phylo) → PhyloNet
Quick utility to convert a Phylo object to a network.
"""
function phylo_to_net(tree)
    return(PhyloNetworks.readnewick(Phylo.outputtree(tree,Newick())))
end

"""
    net_to_phylo(net:PhyloNetwork) → Phylo
Quick utility to convert a PhyloNetwork object to a Phylo tree.
"""
function net_to_phylo(net)
    return(Phylo.parsenewick(PhyloNetworks.writenewick(net)))
end


"""
    rescale_tree(tree:Phylo,scale) -> Phylo
All branch lengths in `tree` are rescaled by multiplying by scale.
"""

function rescale_tree!(tree,scale)
    for e in getbranches(tree)
        e.length = e.length * scale
    end
end

"""
    tip_overlap(tree1::Phylo,tree2::Phylo) -> [String]
returns terminal nodes present in both trees
"""
function tip_overlap(tree1,tree2)
    return(intersect(getleafnames(tree1),getleafnames(tree2)))
end