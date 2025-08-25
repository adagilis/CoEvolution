


function calculate_df(gene_id,vcf_loc,pop_file,tree_file)
    vcf=data_dir*"vcfs/"*gene_id*".vcf.gz"
    isfile(vcf) || make_vcf(gene_id,vcf_loc)
    out = data_dir*"introgression/"*gene_id
    run(`Dsuite Dtrios $vcf $pop_file -t $tree_file -o $out`)
end

function make_vcf(gene_id,vcf_loc)
    id = findfirst(gene_table.gene .== gene_id)
    region = "chr"*gene_table.chr[id]*":"*string(gene_table.pos_1[id])*"-"*string(gene_table.pos_2[id])
    out = data_dir*"vcfs/"*gene_id*".vcf.gz"
    run(`bcftools view -r $region -Oz -o $out $vcf_loc`)
end

function combine_d_stats()
    int_files = filter(endswith("_tree.txt"),readdir(data_dir*"introgression/",join=true))
    genes = replace.(filter(endswith("_tree.txt"),readdir(data_dir*"introgression/",join=false)),r"_tree.txt"=>s"")
    
end