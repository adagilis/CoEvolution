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
    df = DataFrame(:from_x=>from_x,:from_y=>from_y,:to_x=>to_x,:to_y=>to_y)
    network_plot = GLMakie.linesegments(map(floc,eachrow(df)))
    #plot nodes
    GLMakie.scatter!(df[:,1],df[:,2])
    return(network_plot)
end

function floc(i);((i[1],i[3]),(i[2],i[4]));end