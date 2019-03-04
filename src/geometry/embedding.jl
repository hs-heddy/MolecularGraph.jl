#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    graphdistembedding


"""
    graphdistembedding(graph::UDGraph) -> Cartesian2D

Compute 2D embedding based on the graph distance.
"""
function graphdistembedding(graph::UDGraph)
    n = nodecount(graph)
    D = distancematrix(graph)
    H = LinearAlgebra.I - fill(1 / n, (n, n))
    G = - 0.5 * H * D * H # Gram matrix
    F = LinearAlgebra.eigen(G)
    od = sortperm(F.values, rev=true)
    xidx = findall(i -> i == 1, od)[1]
    yidx = findall(i -> i == 2, od)[1]
    coords = zeros(n, 2)
    coords[:, 1] = F.vectors[:, xidx] * sqrt(abs(F.values[xidx]))
    coords[:, 2] = F.vectors[:, yidx] * sqrt(abs(F.values[yidx]))
    return Cartesian2D(coords)
end
