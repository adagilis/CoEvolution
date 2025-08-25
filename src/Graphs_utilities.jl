using Graphs
using SimpleWeightedGraphs
using GraphCommunities
using SGtSNEpi

function construct_ERC_graph(ERC)
    nodes=sort(unique(union(ERC.i,ERC.j)))
    convert=Dict(nodes .=> 1:length(nodes))
    sources = [convert[x] for x in ERC.i]
    destinations = [convert[x] for x in ERC.j]
    weights = abs.(ERC.r)
    g = SimpleWeightedGraph(sources,destinations,weights)
    return(g)
end

function Louvain_cluster(graph)
    communities = compute(Louvain(),SimpleGraph(graph))
    return(communities)
end
