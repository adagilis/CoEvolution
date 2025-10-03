using BioMart

import BioMart: Dataset, Filters, Attributes


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


function num_rbh(gene)
    length(getleafnames(read_tree(data_dir*"trees/"*gene*".treefile")))
end