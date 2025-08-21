function chr_order_pos(position,chr,chr_order,chr_lengths)
    for c in 2:length(chr_order)
        add_dist=sum(chr_lengths[1:(c-1)])
        current_c=chr_order[c]
        id = findall(chr .== current_c)
        position[id] = position[id] .+ add_dist
    end
    return(position)
end


function graph_plot(graph,dataframe)
    
    #plot nodes
    cols=scale_colors(dataframe.community)
    sizes=scale_sizes(dataframe.degree)
    network_plot= @df dataframe scatter(:x,:y,
    c=cols,
    legend=false,
    size=sizes;
    axis=([],false),
    markerstrokewidth=0)

    #plot edges
    id_from = missing
    id_to = missing
    df.from_x = dataframe.x[id_from]
    df.from_y = dataframe.y[id_from]
    df.to_x = dataframe.x[id_to]
    df.to_y = dataframe.y[id_to]
    plot!(
        [([i.from_x,i.to_x],[i.from_y,i.to_y])  for i in eachrow(df)],
        c=:black,alpha=0.5,legend=false)
end