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

function calculate_ERC(gene_id_1::Int64,gene_id_2::Int64,species_tree::Phylo,cutoff=-5)
    scale_rates=deepcopy(species_tree)
    gene1 = open(parsenewick,trees[gene_id_1])
    gene2 = open(parsenewick,trees[gene_id_2])
    tips_both=intersect(getleafnames(tree_i),getleafnames(tree_j))
    keeptips!(tree_i,tips_both)
    keeptips!(tree_j,tips_both)
    keeptips!(tree_tot,tips_both)
    li = [e.length for e in getbranches(tree_i)]
    lj = [e.length for e in getbranches(tree_j)]
    lt = [e.length for e in getbranches(tree_tot)]
    li = li ./ lt
    lj = lj ./ lt
    kept_edges = intersect(findall(<(cutoff),li),findall(<(cutoff),lj))
    if(length(kept_edges)>4)
        li = (li .- mean(li)) ./ (std(li))
        lj = (lj .- mean(lj)) ./ (std(lj))
        r2=cor(li[kept_edges],lj[kept_edges])
        pval=pvalue(CorrelationTest(li[kept_edges],lj[kept_edges]))
    else 
        r2=0
        pval=1
    end
    return(Dict(:i => trees[i],:j => trees[j],:n_edges=>length(kept_edges),:r2 => r2,:pval => pval))
end

#TODO - parallelize

function runERC()
    local ERC_res=[]
        for i in range(1,length(trees)-1)
            for j in range(i+1,length(trees))
                push!(ERC_res,ERC(i,j,total,5))
            end
        end
    return(ERC_res)
end
