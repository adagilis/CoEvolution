function chr_order_pos(position,chr,chr_order,chr_lengths)
    for c in 2:length(chr_order)
        add_dist=sum(chr_lengths[1:(c-1)])
        current_c=chr_order[c]
        id = findall(chr .== current_c)
        position[id] = position[id] .+ add_dist
    end
    return(position)
end
