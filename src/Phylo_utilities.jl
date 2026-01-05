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
function rescale_tree(tree,scale)
    tree_c = deepcopy(tree)
    for e in getbranches(tree_c)
        e.length = e.length * scale
    end
    return(tree_c)
end


"""
    rescale_tree!(tree:Phylo,scale)
All branch lengths in `tree` are rescaled by multiplying by scale, modifies existing tree
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
Takes a tree and a node and returns the name of all terminal descendents of that node
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

"""
    run_astral(trees) -> ASTRAL4 consensus tree
"""
function run_astral(trees;name="astral_consensus.newick")
    treefile=data_dir*"trees/concat_trees.tre"
    for t in trees
        run(pipeline(`cat $t`),stdout=treefile,append=true)
    end
    run(`astral4 -i $treefile -o $data_dir/trees/$name --root $outgroup`)
    run(`rm $treefile`)
end

"""
    run_iqtree(seq) -> iqtree gene tree
Takes a sequence name and runs iqtree on an `aligned.fasta` file in the `aligned` folder.
"""
function run_iqtree(seq)
    cmd = `iqtree2 -s $seq -ntmax 4 -quiet`
    run(pipeline(cmd;stderr=devnull))
    treefile_old = seq*".treefile"
    treefile_new = replace(replace(seq,r".aligned.fasta"=>s".treefile"),r"/aligned/"=>s"/trees/")
    cmd2 = `mv $treefile_old $treefile_new`
    run(cmd2)
    #and clean up
    try
        run(`rm $seq.mldist`)
        run(`rm $seq.model.gz`)
        run(`rm $seq.ckp.gz`)
        run(`rm $seq.iqtree`)
        run(`rm $seq.bionj`)
        run(`rm $seq.log`)
    catch
    end
end

"""
    run_constrained_iqtree(seq,species_tree)
Runs iqtree constraining tree topology to species tree. Useful to compare to ERC2.0.
"""
function run_constrained_iqtree(seq,species_tree)
    taxa=taxa_in_alignment(seq)
    cnstrn=deepcopy(species_tree)
    keeptips!(cnstrn,taxa)
    writenewick(seq*".constraint",cnstrn)
    cmd = `iqtree2 -s $seq -ntmax 4 -quiet -t $seq.constraint`
    run(pipeline(cmd;stderr=devnull))
    treefile_old = seq*".treefile"
    treefile_new = replace(replace(seq,r".aligned.fasta"=>s".treefile"),r"/aligned/"=>s"/trees/")
    cmd2 = `mv $treefile_old $treefile_new`
    run(cmd2)
    #and clean up
    try
        run(`rm $seq.mldist`)
        run(`rm $seq.model.gz`)
        run(`rm $seq.ckp.gz`)
        run(`rm $seq.iqtree`)
        run(`rm $seq.bionj`)
        run(`rm $seq.log`)
        run(`rm $seq.constraint`)
    catch
    end
end


"""
    taxa_in_alignment(file)
Returns the taxa that exists in an aligned fasta file. 
"""
function taxa_in_alignment(file)
    alignment=open(file,"r")
    out=[]
    while(!eof(alignment))
        l=readline(alignment)
        startswith(l,">") && append!(out,[l[2:length[l]]])
    end
    return(out)
end