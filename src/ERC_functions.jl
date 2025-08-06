using Phylo
using FLoops
using DataFrames
using Statistics
using HypothesisTests
using CSV
using ProgressMeter
using Arrow

include("Phylo_utilities.jl")
include("Stats_utilities.jl") 
"""
    calculate_ERC(gene1,gene2,species_tree,cutoff) → [cor,p-value]
Function which takes two gene trees and a species tree, and returns the Evolutionary Rate Correlation between the two genes.
ERC is defined as the correlation of Z-scores in evolutionary rates for both trees.

"""

function calculate_ERC(id1,id2,tree1,tree2,species_tree::RootedTree;cutoff=1)
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
        all_branches = innerjoin(innerjoin(scale_rates,branches_1,on=:node),branches_2,on=:node)
        dropmissing!(all_branches)
        rates_1 = all_branches[:,"bl_1"] ./ all_branches[:,"bl_sp"]
        rates_2 = all_branches[:,"bl_2"] ./ all_branches[:,"bl_sp"]
        z_1 = zscores(rates_1)
        z_2 = zscores(rates_2)
        kept_edges = intersect(findall(<(cutoff),abs.(z_1)),findall(<(cutoff),abs.(z_2)))
        if(length(kept_edges)>3)
            #z_1 = zscores(rates_1[kept_edges])
            #z_2 = zscores(rates_2[kept_edges])
            r = cor(z_1[kept_edges],z_2[kept_edges])
            pval = pvalue(CorrelationTest(z_1[kept_edges],z_2[kept_edges]))
        else 
            #Too few edges with values below cutoff
            r = missing
            pval = missing
        end
    else
        #No tip overlap
        r = missing
        pval = missing
        kept_edges = []
    end
    return(Dict(:i => id1,:j => id2,:n_edges=>length(kept_edges),:r => r,:pval => pval))
end



"""
    runERC_files(trees,species_tree) → DataFrame{[i,j,branches,cor,p-value]}
Calculate a set of ERC values given a list of gene trees `trees` and a species tree `species_tree`. Returns DataFrame object with five columns. 
"""

function runERC_files(trees,species_tree;cutoff=5,ckp_freq=1000000,ckp_file=nothing)
    num_comp = binomial(length(trees),2)
    total = collect(1:num_comp)
    local ERC_res=DataFrame(missings(Float64,num_comp,5),[:i,:j,:n_edges,:r,:pval])
    ERC_res.i = reduce(vcat,[repeat([x],inner=length(trees)-x) for x in 1:length(trees)])
    ERC_res.j = reduce(vcat,[collect(x:12516) for x in 2:12516])
    done = 0
    if (!isnothing(ckp_file))
        if isfile(ckp_file)
            ERC_res = DataFrame(Arrow.Table(ckp_file))
            total = setdiff(total,findall(completecases(ERC_res)))
            done = length(findall(completecases(ERC_res)))
        end
    end
    p = Progress(num_comp,desc="Calculating ERC scores:")
    @floop ThreadedEx() for index in total
        i = ERC_res.i[index]
        j = ERC_res.j[index]
        ti = read_tree(trees[i])
        tj = read_tree(trees[j])
        try
            ERC_res[index,:] = calculate_ERC(i,j,ti,tj,species_tree;cutoff=cutoff)
        catch
            #println("ERC failed to calculate for $(i), $(j)")
            #If this happens - something went wrong! We keep r2 different from 0 to be able to quantify how frequently
            ERC_res[index,:] = Dict(:i => i,:j => j,:n_edges =>missing,:r => missing,:pval =>missing)
        end
        done += 1
        if (done % ckp_freq == 0 & !isnothing(ckp_file))
            Arrow.write(ckp_file,ERC_res)
        end
        next!(p)
    end
    return(ERC_res)
end


"""
    runERC_collection(trees,species_tree::Phylo;cutoff=5) -> DataTable[i,j,edges_kept,r2,pval]
    calculates ERC values using the `calculate_ERC` function for all pairs of trees in the trees object.
"""

function runERC_collection(trees,species_tree;cutoff=5)
    num_comp = binomial(length(trees),2)
    p = Progress(num_comp,desc="Calculating ERC scores:")
    local ERC_res=DataFrame(zeros(num_comp,5),[:i,:j,:n_edges,:r,:pval])
    @floop ThreadedEx() for (i,j) in Iterators.product(1:(length(trees)-1),1:length(trees))
        if(j>i)
            index = index_func(i,j,length(trees))
            try
                ERC_res[index,:] = calculate_ERC(i,j,trees[i],trees[j],species_tree;cutoff=cutoff)
            catch
                #println("ERC failed to calculate for $(i), $(j)")
                ERC_res[index,:] = Dict(:i => i,:j => j,:n_edges =>0,:r => -1,:pval =>-1)
            end
            next!(p)
        end
    end
    return(ERC_res)
end