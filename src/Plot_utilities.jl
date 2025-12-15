function chr_order_pos(position,chr,chr_order,chr_lengths)
    for c in 2:length(chr_order)
        add_dist=sum(chr_lengths[1:(c-1)])
        current_c=chr_order[c]
        id = findall(chr .== current_c)
        position[id] = position[id] .+ add_dist
    end
    return(position)
end


function graph_plot(graph,dataframe;cols=:blue,sizes=1,linealpha=0.05)
    #plot edges
    A = adjacency_matrix(graph)
    nz = findnz(triu(A))
    id_from = nz[1]
    id_to = nz[2]
    weights = nz[3]
    from_x = dataframe[:,1][id_from]
    from_y = dataframe[:,2][id_from]
    to_x = dataframe[:,1][id_to]
    to_y = dataframe[:,2][id_to]
    df = DataFrame(:from_x=>from_x,:from_y=>from_y,:to_x=>to_x,:to_y=>to_y,:weight=>weights)
    network_plot = GLMakie.linesegments(map(floc,eachrow(df)))
    #plot nodes
    GLMakie.scatter!(df[:,1],df[:,2])
    return(network_plot)
end

"""
    sig2mat
Function which takes an ERC result data-table, and returns an n x n matrix, where each entry (i,j) is 1-ERC^2 value for genes i and j
"""
function sig2mat(ERC;default=0)
    scores = 1.0 .- ERC.r.^2
    genes = union(ERC.i,ERC.j)
    retmat=fill(default,length(genes),length(genes))
    for g in 1:length(genes)
        id_g_i = findall(ERC.i .==g)
        id_g_j = findall(ERC.j .==g)
        id_2_i = indexin(ERC.j[id_g_i],genes)
        id_2_j = indexin(ERC.i[id_g_j],genes)
        retmat[id_g_i,id_2_i] = scores[id_g_i]
        retmat[id_g_j,id_2_j] = scores[id_g_j]
    end
    return(retmat)
end