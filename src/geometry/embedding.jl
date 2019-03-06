#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    graphdistembedding


"""
    graphdistembedding(graph::UDGraph) -> Coordinates

Compute cartesian embedding based on the graph distance and returns
`Cartesian2D` or `Cartesian3D` which is more appropriate. If `force3d` is true,
it always returns `Cartesian3D`.
"""
function graphdistembedding(graph::UDGraph; yzfactor=2.0, force3d=false)
    n = nodecount(graph)
    D = distancematrix(graph)
    H = LinearAlgebra.I - fill(1 / n, (n, n))
    G = - 0.5 * H * D * H # Gram matrix
    F = LinearAlgebra.eigen(G)
    display(F)
    od = sortperm(real.(F.values), rev=true)
    xidx = od[1]
    yidx = od[2]
    zidx = od[3]
    coords = zeros(n, 3)
    xv = abs(F.values[xidx])
    yv = abs(F.values[yidx])
    zv = abs(F.values[zidx])
    println([xv, yv, zv])
    coords[:, 1] = real.(F.vectors[:, xidx]) * sqrt(xv)
    coords[:, 2] = real.(F.vectors[:, yidx]) * sqrt(yv)
    coords[:, 3] = real.(F.vectors[:, zidx]) * sqrt(zv)
    if !force3d && yv / zv > yzfactor
        return Cartesian2D(coords[:, 1:2])
    else
        return Cartesian3D(coords)
    end
end
