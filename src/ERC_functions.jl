using Phylo
using FLoops
using DataFrames
using Statistics
using HypothesisTests
using CSV
using ProgressMeter

include("Phylo_utilities.jl") 
"""
    calculate_ERC(gene1,gene2,species_tree,cutoff) → [cor,p-value]
Function which takes two gene trees and a species tree, and returns the Evolutionary Rate Correlation between the two genes.
ERC is defined as the correlation of Z-scores in evolutionary rates for both trees.

"""

function calculate_ERC(id1,id2,tree1,tree2,species_tree::RootedTree;cutoff=5)
    template = deepcopy(species_tree)
    gene1 = deepcopy(tree1)
    gene2 = deepcopy(tree2)
    tips_both=tip_overlap(gene1,gene2)
    if(length(tips_both) >= 3)
        keeptips!(template,tips_both)
        keeptips!(gene1,tips_both)
        keeptips!(gene2,tips_both)        
        scale_rates = breakdown_tree(template)
        rename!(scale_rates,:bl => :bl_sp)
        branches_1 = breakdown_tree(gene1)
        rename!(branches_1,:bl => :bl_1)
        branches_2 = breakdown_tree(gene2)
        rename!(branches_2,:bl => :bl_2)
        all_branches = leftjoin(leftjoin(scale_rates,branches_1,on=:node),branches_2,on=:node)
        dropmissing!(all_branches)
        rates_1 = all_branches[:,"bl_1"] ./ all_branches[:,"bl_sp"]
        rates_2 = all_branches[:,"bl_2"] ./ all_branches[:,"bl_sp"]
        kept_edges = intersect(findall(<(cutoff),rates_1),findall(<(cutoff),rates_2))
        if(length(kept_edges)>3)
            z_1 = zscore(rates_1)
            z_2 = zscore(rates_2)
            r2 = cor(z_1[kept_edges],z_2[kept_edges])
            pval = pvalue(CorrelationTest(z_1[kept_edges],z_2[kept_edges]))
        else 
            #Too few edges with values below cutoff
            r2 = 0.0
            pval = 1.0
        end
    else
        #No tip overlap
        r2 = 0.0
        pval = 1.0
        kept_edges = []
    end
    return(Dict(:i => id1,:j => id2,:n_edges=>length(kept_edges),:r2 => r2,:pval => pval))
end



"""
    runERC_files(trees,species_tree) → DataFrame{[i,j,branches,cor,p-value]}
Calculate a set of ERC values given a list of gene trees `trees` and a species tree `species_tree`. Returns DataFrame object with five columns. 
"""

function runERC_files(trees,species_tree;cutoff=5)
    num_comp = binomial(length(trees),2)
    p = Progress(num_comp,desc="Calculating ERC scores:")
    local ERC_res=DataFrame(zeros(num_comp,5),[:i,:j,:n_edges,:r2,:pval])
    @floop ThreadedEx() for (i,j) in Iterators.product(1:(length(trees)-1),1:length(trees))
        if(j>i)
            ti = read_tree(trees[i])
            tj = read_tree(trees[j])
            index = index_func(i,j,length(trees))
            try
                ERC_res[index,:] = calculate_ERC(i,j,ti,tj,species_tree;cutoff=5)
            catch
                #println("ERC failed to calculate for $(i), $(j)")
                #If this happens - something went wrong! We keep r2 different from 0 to be able to quantify how frequently
                ERC_res[index,:] = Dict(:i => i,:j => j,:n_edges =>0,:r2 => -1,:pval =>-1)
            end
            next!(p)
        end
    end
    return(ERC_res)
end

function runERC_collection(trees,species_tree;cutoff=5)
    num_comp = binomial(length(trees),2)
    p = Progress(num_comp,desc="Calculating ERC scores:")
    local ERC_res=DataFrame(zeros(num_comp,5),[:i,:j,:n_edges,:r2,:pval])
    @floop ThreadedEx() for (i,j) in Iterators.product(1:(length(trees)-1),1:length(trees))
        if(j>i)
            index = index_func(i,j,length(trees))
            try
                ERC_res[index,:] = calculate_ERC(i,j,trees[i],trees[j],species_tree;cutoff=5)
            catch
                #println("ERC failed to calculate for $(i), $(j)")
                ERC_res[index,:] = Dict(:i => i,:j => j,:n_edges =>0,:r2 => -1,:pval =>-1)
            end
            next!(p)
        end
    end
    return(ERC_res)
end

function index_func(i,j,n)
    k = 0
    if(i>1)
        k = (i-1)*n-sum(1:i .- 1)
    end
    return(k+j-i)
end

function zscore(x)
    return((x .- mean(x))./ std(x))
end

function get_branch_lengths(tree)
    return([e.length for e in getbranches(tree)])
end

function reorder_tree(tree::RootedTree,tips)
    sort!(keeptips!(tree,tips))
    return(sort!(Phylo.parsenewick(Phylo.outputtree(tree,Newick()))))
end

function cladewise_order!(tree)
    local tmp = phylo_to_net(tree)
    PhyloNetworks.cladewiseorder!(tmp)
    return(net_to_phylo(tmp))
end

#Quick breakdown of branches so that each branch is designated by the node it terminates in.
"""
    breakdown_tree(tree::Phylo) → DataFrame(bl::Float64,node::String)
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