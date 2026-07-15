using BioMart
using CSV
using DataFrames
import BioMart: Dataset, Filters, Attributes
using MultipleTesting: adjust, BenjaminiHochberg


"""
    find_name(gene,dataset) -> BioMart annotations for a sequence.
"""
function find_name(gene,dataset)
    res = BioMart.query(
        Dataset(dataset),
        Filters(refseq_peptide=gene),
        Attributes(
            "refseq_peptide",
            "ensembl_gene_id",
            "external_gene_name",
            "chromosome_name",
            "start_position",
            "end_position",
            "description"
        )
    )
    rename!(res,:var"RefSeq peptide ID"=>:gene,
                :var"Gene stable ID"=>:stable_id,
                :var"Gene name"=>:name,
                :var"Chromosome/scaffold name"=>:chr,
                :var"Gene start (bp)"=>:pos_1,
                :var"Gene end (bp)"=>:pos_2,
                :var"Gene description"=>:description)
    return(res)
end

"""
    num_rbh(gene) -> number of taxon in tree for a gene
"""
function num_rbh(gene)
    length(getleafnames(read_tree(data_dir*focal*"/trees/"*gene*".treefile")))
end

"""
    open_GO(file) -> DataFrame
    Function to read in .gfa files with appropriate column headers. Tested only on v. of .GFA.
"""
function open_GO(file)
    table = CSV.read(file,DataFrame,header=false,comment="!")
    rename!(table,:Column1=>:DB,
    :Column2=>:DB_object_id,
    :Column3=>:DB_object_symbol,
    :Column4=>:Qualifier,
    :Column5=>:GO_ID,
    :Column6=>:DB_reference,
    :Column7=>:Evidence_Code,
    :Column8=>:WithOrFrom,
    :Column9=>:Aspect,
    :Column10=>:DB_object_name,
    :Column11=>:DB_object_synonym,
    :Column12=>:DB_object_type,
    :Column13=>:Taxon,
    :Column14=>:Date,
    :Column15=>:Assigned_by,
    :Column16=>:Annotation_Extension,
    :Column17=>:Gene_Product_Form_ID)
    return(table)
end

"""
    gene_GO(stable_id;score=1) -> DataFrame(:GO,:score)
Given a GO_table object exists, return a list of unique GO IDs for a gene. Returns a DataFrame with scores for each GO ID. Helper function to summarize GO term propogation through ERC network.
"""
function gene_GO(stable_id,GO_table;score=1)
    if !ismissing(stable_id)
        subGO = filter(:DB_object_id=> x -> x==stable_id,GO_table)
        uniGO = unique(subGO.GO_ID)
        return(DataFrame(:GO=>uniGO,:score=>score))
    else
        return(DataFrame(:GO=>missing,:score=>missing))
    end
end


"""
    ERC_GO_extend(gene_id,ERC,GO_db,gene_table;back_GO=missing) -> DataFrame(:GO_term,:score,:occurence,:adjusted_score)
For a given gene with internal gene_id and set of evolutionary rate correlations (ERC - either full list or only significant), returns a list of GO terms of interactions partners, weighted by their co-evolutionary score. 
Recommended to havea pre-calculated backround GO expectation (`go_null` function), as well as a table defining the GO terms for each gene (`GO_db`).
"""
function ERC_GO_extend(gene_id,ERC,gene_table,GO_db,back_GO)
    back_dict = Dict(back_GO.GO_ID .=> 1:length(back_GO.GO_ID))
    subERC = filter([:i,:j] => (i,j) -> i==gene_id || j==gene_id,ERC)
    partners = setdiff(unique(hcat(subERC.i,subERC.j)),[gene_id])
    idx = indexin(partners,gene_table.gene)
    partners = partners[findall((!isnothing).(idx))]
    idx = idx[findall((!isnothing).(idx))]
    fbid = gene_table.flybase[idx]
    go_table = reduce(vcat,[gene_GO(fbid[x],GO_db;score=subERC.fERC[findfirst(subERC.i .== partners[x] .|| subERC.j .== partners[x])]) for x in 1:length(fbid)])
    go_table = disallowmissing(go_table[completecases(go_table),:])
    if size(go_table) != (0,0)
        df = groupby(go_table,:GO)
        uniGO = unique(go_table.GO)
        MUtests = [length(df[df.keymap[(go,)]].score) >2 && MannWhitneyUTest(df[df.keymap[(go,)]].score,back_GO.expected[back_dict[go]]) for go in uniGO]
        ids = findall((!isa).(MUtests,Bool))
        ret = DataFrame(:GO => uniGO[ids],
            :pval=>pvalue.(MUtests[ids]),
            :exp_fERC=>mean.(back_GO.expected[[back_dict[g] for g in uniGO[ids]]]),
            :obs_fERC=>combine(df,:score=>mean).score_mean[ids],
            :pval_BH=>adjust(pvalue.(MUtests[ids]),BenjaminiHochberg()))
        return(sort(ret,:pval))
    else
        return(DataFrame(:GO=>missing,:pval=>missing,:exp=>0,:obs_mean=>0,:pval_BH=>missing))
    end
end

"""
    go_null(GO_table,gene_table) -> DataFrame(:GO_term,:expected_score)
Generates a set of GO terms for an ERC network, with weights relative to the average ERC value for all genes with that GO term. Useful as a baseline for GO term propogation, used for t-tests later.
"""
function go_null(GO_table,gene_table)
    uniGO = unique(GO_table.GO_ID)
    mean_scores = [collect(skipmissing(filter(:flybase=> g ->!ismissing(g) && g ∈  GO_table.DB_object_id[GO_table.GO_ID .== GO],gene_table).mean_fERC)) for GO in uniGO]
    return(DataFrame(:GO_ID=>uniGO,:expected=>mean_scores))
end

"""
    go2gene(go,GO_table,gene_table)
Return the gene ids for all genes in any GO category. Useful to plot GO categories across the network.
"""
function go2gene(go,GO_table,gene_table)
    fbids = filter(:GO=>g->g==go,GO_table).DB_object_id
    geneids = filter(:flybase=>fb-> fb ∈ fbids,gene_table).gene
    return(geneids)
end

"""
    genes_w_GO(GO_ID,GO_table,gene_table) -> gene_ids
function to take a GO ID and return a set of gene ids in your gene_table corresponding to that GO term. Useful for visualization.
"""
function genes_w_GO(GO_ID,GO_table,gene_table)
    fbids = unique(filter(:GO_ID=>gid -> gid==GO_ID,GO_table).DB_object_id)
    gene_ids = filter(:flybase=>fid -> !ismissing(fid) && fid ∈ fbids,gene_table).gene
    return(gene_ids)
end

"""
    annotate_GO_table!(back_GO,go.obo) => back_GO_w_annotation
    Adds name and namespace to each go term in DB, so that results are human readable/interpretable.
"""
function annotate_GO_table!(back_GO,go_obo)
    ret = [fetch_name_GO(go,go_obo) for go in back_GO.GO_ID]
    names, types = map(x->getindex.(ret,x),1:2)
    back_GO.type .= types
    back_GO.desc .= names
end

"""
fetch_name_GO(go_term,go_obo)
Given the location of an OBO file, returns what the description of a go term is and which namespace (BP,MF,CC) it's in.
"""
function fetch_name_GO(go_term,go_obo)
    cmd = `grep "id: $go_term" -A 2 $go_obo`
    ret=split(read(cmd,String),"\n")[2:3]
    name = replace(ret[1],"name: "=>"")
    type = replace(ret[2],"namespace: "=>"")
    return((name,type))
end
