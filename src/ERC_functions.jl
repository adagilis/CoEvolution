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
                #Too few edges with values below cutoff
                r2 = 0.0
                pval = 1.0
            end
        else
            #Insufficient edge overlap
            r2 = 0.0
            pval = 1.0
            kept_edges=[]
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
    return(parsenewick(Phylo.outputtree(tree,Newick())))
end