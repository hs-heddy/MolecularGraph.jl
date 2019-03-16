#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    merge


struct Cartesian2DMerger
    groups::Vector{Cartesian2DSubset}
    links::Dict{Int,Dict{Int,Set{Int}}}

    Cartesian2DMerger() = new([], Dict())
end


"""
    merge(merger::Cartesian2DMerger) -> Cartesian2D

Merge 2D cartesian coordinates subsets.
"""
function merge!(merger::Cartesian2DMerger)
    merger.links

end


function dfs!(merger::Cartesian2DMerger, n::Int)
    p = merger.pred[n]

    link = merger.links[p][n]
    if length(link) == 1
        # point merge
        translate =
        rotation =
    elseif length(link) == 2
        # line merge
        translate =
        rotation = 
    end
    for nbr in neighborkeys(merger.groups, n)
        merger.pred[n] == nbr && continue
        dfs!(merger, nbr)
    end
end
