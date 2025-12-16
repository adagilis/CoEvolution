using BioMart
using CSV
using DataFrames
import BioMart: Dataset, Filters, Attributes


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
            "flybase_gene_id",
            "chromosome_name",
            "start_position",
            "end_position",
            "description"
        )
    )
    rename!(res,:var"RefSeq peptide ID"=>:gene,
                :var"Gene stable ID"=>:stable_id,
                :var"Gene name"=>:name,
                :var"FlyBase gene ID"=>:flybase,
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
    length(getleafnames(read_tree(data_dir*"trees/"*gene*".treefile")))
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
function gene_GO(stable_id;score=1)
    if !ismissing(stable_id)
        subGO = filter(:DB_object_id=> x -> x==stable_id,GO_table)
        uniGO = unique(subGO.GO_ID)
        return(DataFrame(:GO=>uniGO,:score=>score))
    else
        return(DataFrame(:GO=>missing,:score=>missing))
    end
end


"""
    ERC_GO_extend(gene_id,ERC) -> DataFrame(:GO_term,:score,:occurence,:adjusted_score)
    Requires `gene_table` and `GO_table` to exist. For a given gene and set of evolutionary rate correlations, returns a list of GO terms of interactions partners, weighted by their co-evolutionary score. 
    Definitely needs a better approach to capture null expectation. Simplest approach would be to perform GO-enrichment among significant terms, but the current approach allows weighting by the relative ERC score, which is harder to capture with GO-enrichment.
"""
function ERC_GO_extend(gene_id,ERC)
    subERC = filter([:i,:j] => (i,j) -> i==gene_id || j==gene_id,ERC)
    partners = setdiff(unique(hcat(subERC.i,subERC.j)),[gene_id])
    fbid = gene_table.flybase[indexin(partners,gene_table.gene)]
    go_table = reduce(vcat,[gene_GO(fbid[x];score=subERC.fERC[findfirst(subERC.i .== partners[x] .|| subERC.j .== partners[x])]) for x in 1:length(fbid)])
    if size(go_table) != (0,0)
        df = combine(groupby(go_table,:GO),:score=> sum,nrow=>:occurence)
        df = df[completecases(df),:]
        df.adjust= df.score_sum ./ df.occurence
        return(sort(df,:adjust))
    else
        return(DataFrame(:GO=>missing,:score_sum=>missing,:occurence=>0,:adjust=>0))
    end
end

"""
    go_backbround(ERC) -> DataFrame(:GO_term,:score,:occurence)
    Generates a set of GO terms for an ERC network, with weights relative to the average ERC value for all genes with that GO term. Useful as a baseline for GO term propogation
"""
function go_background(ERC)
    genes=unique(hcat(ERC.i,ERC.j))
    ave_ERC = [mean(filter([:i,:j] => (i,j) -> i == g || j == g,ERC).fERC) for g in genes]
    go_table = reduce(vcat,[gene_GO(genes[x];score=ave_ERC[x]) for x in 1:length(genes)])
    if size(go_table) != (0,0)
        df = combine(groupby(go_table,:GO),:score=> sum,nrow=>:occurence)
        df = df[completecases(df),:]
        df.adjust= df.score_sum ./ df.occurence
        return(sort(df,:adjust))
    else
        return(DataFrame(:GO=>missing,:score_sum=>missing,:occurence=>0,:adjust=>0))
    end
end