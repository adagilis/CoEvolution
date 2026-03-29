using Phylo
using FLoops
using DataFrames
using Statistics
using HypothesisTests
using CSV
using ProgressMeter
using Term
using Term.Tables
using Arrow
import Term: tprint

include("Phylo_utilities.jl")
include("Stats_utilities.jl") 
"""
    calculate_ERC(gene1,gene2,species_tree;cutoff=5) → [cor,p-value]
Function which takes two gene trees and a species tree, and returns the Evolutionary Rate Correlation between the two genes.
ERC is defined as the correlation of relativized evolutionary rates between two trees. We use the approach as defined by Clark, with the default cutoff of 5 chosen based on Steenwyk et al.

Important notes: I think this approach has more to give. There are definitely issues with discounting branches above an arbitrary rate threshold - rapidly evolving genes can get dropped entirely.
Similarly, I don't perform the z-transform as done in both prior approaches because it's simply redundant (cov(zx,zy)=cor(x,y)=cov(x,y)=cor(zx,zy)).

In short - there's definitely more to develop this approach, but we're sticking to published methods here for simplicity.
"""
function calculate_ERC(id1,id2,tree1,tree2,species_tree::RootedTree;cutoff=5,min_shared=3,min_edges=4)
    template = deepcopy(species_tree)
    gene1 = deepcopy(tree1)
    gene2 = deepcopy(tree2)
    tips_both=tip_overlap(gene1,gene2)
    shared_tips = length(tips_both)
    if(length(tips_both) >= min_shared)
        keeptips!(template,tips_both)
        keeptips!(gene1,tips_both)
        keeptips!(gene2,tips_both)        
        scale_rates = breakdown_tree(template)
        rename!(scale_rates,:bl => :bl_sp)
        branches_1 = breakdown_tree(gene1)
        branches_1 = branches_1[findall(branches_1.bl .> 0),:]
        rename!(branches_1,:bl => :bl_1)
        branches_2 = breakdown_tree(gene2)
        branches_2 = branches_2[findall(branches_2.bl .> 0),:]
        rename!(branches_2,:bl => :bl_2)
        all_branches = innerjoin(innerjoin(scale_rates,branches_1,on=:node),branches_2,on=:node)
        dropmissing!(all_branches)
        rates_1 = all_branches[:,"bl_1"] ./ all_branches[:,"bl_sp"] 
        rates_2 = all_branches[:,"bl_2"] ./ all_branches[:,"bl_sp"]
        kept_edges = intersect(findall(<(cutoff),rates_1),findall(<(cutoff),rates_2))
        if(length(kept_edges)>min_edges)
            corTest = CorrelationTest(rates_1[kept_edges],rates_2[kept_edges])
            r = corTest.statistic
            fERC = fisher_trans(r,length(kept_edges))
            pval = pvalue(corTest)
        else 
            #Too few edges with values below cutoff
            r = missing
            fERC = missing
            pval = missing
        end
    else
        #No tip overlap
        r = missing
        fERC = missing
        pval = missing
        kept_edges = []
    end
    return(Dict(:i => id1,:j => id2,:n_edges=>length(kept_edges),:r => r,:fERC => fERC,:pval => pval,:shared_tips => shared_tips))
end



"""
    runERC_files(trees,species_tree) → DataFrame{[i,j,branches,cor,p-value]}
Calculate a set of ERC values given a list of gene trees `trees` and a species tree `species_tree`. Returns DataFrame object with five columns. 
"""
function runERC_files(trees,species_tree;cutoff=5,min_shared=3,min_edges=4)
    tprint("{blue}fERC score calculation:{/blue}")
    println("Using "*string(Threads.nthreads())*" threads.")
    num_comp = binomial(length(trees),2)
    println("Will be performing "*@green(string(num_comp))*" comparisons.")
    total = collect(1:num_comp)
    local ERC_res=DataFrame(missings(Float64,num_comp,7),[:i,:j,:n_edges,:r,:fERC,:pval,:shared_tips])
    ERC_res.i = reduce(vcat,[repeat([x],inner=length(trees)-x) for x in 1:length(trees)])
    ERC_res.j = reduce(vcat,[collect(x:length(trees)) for x in 2:length(trees)])
    p = Progress(num_comp,desc="Calculating ERC scores:",showspeed=true)
    @floop ThreadedEx() for index in total
        i = ERC_res.i[index]
        j = ERC_res.j[index]
        try
            ERC_res[index,:] = calculate_ERC(i,j,read_tree(trees[i])read_tree(trees[j]),species_tree;cutoff=cutoff,min_shared=min_shared,min_edges=min_edges)
        catch
            #If this happens - something went wrong! Worth thinking about returning something other than missing to identify these cases
            ERC_res[index,:] = Dict(:i => i,:j => j,:n_edges =>missing,:r => missing,:fERC=>missing,:pval =>missing,:shared_tips=>missing)
        end
        next!(p)
    end
    #completed = length(findall(completecases(ERC_res)))
    #data_table=hcat(["Total Pairs","Completed","Insufficient Data"],[num_comp,completed,num_comp-completed])
    #println(Table(data_table;header=["Class","Number of Pairs"]))
    return(ERC_res)
end


"""
    runERC_collection(trees,species_tree::Phylo;cutoff=5) -> DataTable[i,j,edges_kept,r2,pval]
calculates ERC values using the `calculate_ERC` function for all pairs of trees in the trees object.
"""
function runERC_collection(trees,species_tree;cutoff=5,min_shared=3)
    tprint("{blue}fERC score calculation:{/blue}")
    println("Using "*string(Threads.nthreads())*" threads.")
    num_comp = binomial(length(trees),2)
    println("Will be performing "*@green(string(num_comp))*" comparisons.")
    num_comp = binomial(length(trees),2)
    total = collect(1:num_comp)
    local ERC_res=DataFrame(missings(Float64,num_comp,7),[:i,:j,:n_edges,:r,:fERC,:pval,:shared_tips])
    ERC_res.i = reduce(vcat,[repeat([x],inner=length(trees)-x) for x in 1:length(trees)])
    ERC_res.j = reduce(vcat,[collect(x:length(trees)) for x in 2:length(trees)])
    p = Progress(num_comp,desc="Calculating ERC scores:",showspeed=true)
    @floop ThreadedEx() for index in total
        i = ERC_res.i[index]
        j = ERC_res.j[index]
        try
            ERC_res[index,:] = calculate_ERC(i,j,trees[i],trees[j],species_tree;cutoff=cutoff,min_shared=min_shared)
        catch
            #If this happens - something went wrong! Worth thinking about returning something other than missing to identify these cases
            ERC_res[index,:] = Dict(:i => i,:j => j,:n_edges =>missing,:r => missing,:fERC=>missing,:pval =>missing,:shared_tips=>missing)
        end
        next!(p)
    end
    return(ERC_res) 
end


"""
    function fisher_trans(r,n_edges)-> fisher transformed correlation value
"""
function fisher_trans(r,n_edges)
    return atanh(r) * sqrt(n_edges-3)
end

"""
    function specificity(i,j)-> returns weighted strenght of an interaction.
        Defined as frac{2*fERC_{ij}}{sum_if ERC_{ij} sum_j + fERC_{ij}} 
        Represents how much of the significant fERC signal between these two genes is unique to them.
        Meant to filter out cases where strong interactions with hub genes dominate outputs of analyses - 
        in some cases we may be more interested in less strong signals of evolution that represent more of 
        the co-evolutionary signal for the pair of genes together.

        Note that this could also be calculated using the degree of each gene, without ascribing weights based on fERC.
        That would be, in some ways, more accurately called specificity, but I needed a quick function name.
"""
function specificity(a,b,sigERC)
    sum_a = sum(sigERC.fERC[sigERC.i .==a .|| sigERC.j .==a])
    sum_b = sum(sigERC.fERC[sigERC.i .==b .|| sigERC.j .==b])
    val = sigERC.fERC[sigERC.gid .== join(sort([a,b]),"_")]
    return(2*val/(sum_a+sum_b))
end