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


"""
    get_branch_lengths(tree::Phylo) -> [Float64]
    returns the branch lengths for all branches in `tree`
"""
function get_branch_lengths(tree)
    return([e.length for e in getbranches(tree)])
end

"""
    breakdown_tree(tree::Phylo) → DataFrame(bl::Float64,node::String)
    labels each branch by the node it terminates in, and returns a table with the branch name and lenght
"""
function breakdown_tree(tree)
    sum_table=DataFrame(:bl=>missings(Float64,nbranches(tree)),:node=>missings(String,nbranches(tree)))
    counter=1
    for n in getnodes(tree)
        nodes = getdescendant_leaves(tree,n)
        inbound = getinbound(tree,n)
        if !isnothing(inbound)
            bl = getlength(tree,getinbound(tree,n))
            sum_table[counter,:] = Dict(:bl=>bl,:node=>nodes)
            counter += 1
        end
    end
    return(sum_table)
end

"""
    getdescendant_leaves(tree,n) -> [node_names]
    takes a tree and a node and returns the name of all terminal descendents of that node
"""

function getdescendant_leaves(tree,n)
    #This isn't ideal, but we rely on the fact that internal nodes are relabeled by Phylo.jl to start with "Node"
    all_nodes = Phylo.getdescendants(tree,n)
    if length(all_nodes) > 0
        nodes = filter(!contains("Node"),[i.name for i in all_nodes])
        return(join(map(string,sort!(nodes)),','))
    else
        return(n.name)
    end
end