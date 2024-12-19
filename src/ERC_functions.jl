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
using RCall

function calculate_ERC(gene_id_1::Int64,gene_id_2::Int64,species_tree,cutoff=-5)
    scale_rates=deepcopy(species_tree)
    gene1 = open(parsenewick,trees[gene_id_1])
    gene2 = open(parsenewick,trees[gene_id_2])
    tips_both=intersect(getleafnames(gene1),getleafnames(gene2))
    keeptips!(gene1,tips_both)
    keeptips!(gene2,tips_both)
    keeptips!(scale_rates,tips_both)
    l1 = [e.length for e in getbranches(gene1)]
    l2 = [e.length for e in getbranches(gene2)]
    lt = [e.length for e in getbranches(scale_rates)]
    l1 = l1 ./ lt
    l2 = l2 ./ lt
    kept_edges = intersect(findall(<(cutoff),l1),findall(<(cutoff),l2))
    if(length(kept_edges)>3)
        l1 = (l1 .- mean(l1)) ./ (std(l1))
        l2 = (l2 .- mean(l2)) ./ (std(l2))
        r2=cor(l1[kept_edges],l2[kept_edges])
        pval=pvalue(CorrelationTest(l1[kept_edges],l2[kept_edges]))
    else 
        r2=0
        pval=1
    end
    return(Dict(:i => trees[gene_id_1],:j => trees[gene_id_2],:n_edges=>length(kept_edges),:r2 => r2,:pval => pval))
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
