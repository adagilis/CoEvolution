using Graphs
using SimpleWeightedGraphs
using CommunityDetection

function construct_ERC_graph(ERC_results,trees)
    sources = ERC_results.i
    destinations = ERC_results.j
    weights = ERC_results.r2
    g = SimpleWeightedGraph(sources,destinations,weights)
    return(g)
end