##Make a distance matrix from ERC
function make_dist(ERC)
    keys=genes
    values = [findall(ERC.i .== keys[x] .|| ERC.j .== keys[x]) for x in 1:length(keys)]
    ERC_dict = Dict(keys=>values)
    ret = missings(Float64,[length(keys),length(keys)])
    for (i,j) in Iterators.product(1:length(keys),1:length(keys))
        if i!=j
            out = intersect(ERC_dict[keys[i]],ERC_dict[keys[j]])
            
        end
    end
end