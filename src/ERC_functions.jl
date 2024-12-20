"""
    calculate_ERC(gene1,gene2,species_tree,cutoff) → [cor,p-value]
Function which takes two gene trees and a species tree, and returns the Evolutionary Rate Correlation between the two genes.
ERC is defined as the correlation of Z-scores in evolutionary rates for both trees.

"""

using Phylo
using DataFrames
using Statistics
using HypothesisTests
using CSV

function read_tree(file)
    return(open(parsenewick,file))
end

function calculate_ERC(id1,id2,tree1,tree2,species_tree::RootedTree,cutoff=-5)
    template = deepcopy(species_tree)
    gene1 = deepcopy(tree1)
    gene2 = deepcopy(tree2)
    tips_both=intersect(getleafnames(gene1),getleafnames(gene2))
    if(length(tips_both) >= 3)
        template = reorder_tree(template,tips_both)
        gene1 = reorder_tree(gene1,tips_both)
        gene2 = reorder_tree(gene2,tips_both)
        if(length(gene1.branches)==length(gene2.branches))
            scale_rates = get_branch_lengths(template)
            branches_1 = get_branch_lengths(gene1) ./ scale_rates
            branches_2 = get_branch_lengths(gene2) ./ scale_rates
            kept_edges = intersect(findall(<(cutoff),branches_1),findall(<(cutoff),branches_2))
            if(length(kept_edges)>3)
                z_1 = zscore(branches_1)
                z_2 = zscore(branches_2)
                r2 = cor(z_1[kept_edges],z_2[kept_edges])
                pval = pvalue(CorrelationTest(z_1[kept_edges],z_2[kept_edges]))
            else 
                r2 = 0
                pval = 1
            end
        else
            r2 = 0
            pval = 1
            kept_edges=0
        end
    else
        r2 = 0
        pval = 1
        kept_edges = 0
    end
    return(Dict(:i => id1,:j => id2,:n_edges=>length(kept_edges),:r2 => r2,:pval => pval))
end

#TODO - parallelize

function runERC_files(trees)
    local ERC_res=[]
        for i in range(1,length(trees)-1)
            ti = read_tree(trees[i])
            for j in range(i+1,length(trees))
                tj = read_tree(trees[j])
                push!(ERC_res,calculate_ERC(i,j,ti,tj,species_tree,5))
            end
        end
    return(DataFrame(ERC_res))
end

function runERC_collection(trees)
    local ERC_res=[]
    for i in range(1,length(trees)-1)
        for j in range(i+1,length(trees))
            push!(ERC_res,calculate_ERC(i,j,trees[i],trees[j],species_tree,5))
        end
    end
return(DataFrame(ERC_res))
end

function zscore(x)
    return((x .- mean(x))./ std(x))
end

function get_branch_lengths(tree)
    return([e.length for e in getbranches(tree)])
end

function reorder_tree(tree::RootedTree,tips)
    sort!(keeptips!(tree,tips))
    return(parsenewick(Phylo.outputtree(tree,Newick())))
end