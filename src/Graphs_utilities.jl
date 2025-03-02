using Graphs
using SimpleWeightedGraphs
using GraphCommunities

function construct_ERC_graph(ERC_results)
    sources = Int.(ERC_results.i)
    destinations = Int.(ERC_results.j)
    weights = abs.(ERC_results.r)
    g = SimpleWeightedGraph(sources,destinations,weights)
    return(g)
end

function Louvain_cluster(graph)
    communities = compute(Louvain(),SimpleGraph(graph))
    return(communities)
end
