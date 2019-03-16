#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Embed2DConstraints,
    addgroupconstraint!


struct Embed2DConstraints
    group::Vector{PointSet}
    distance::Vector{Tuple{Int,Int,Float64}}
    angle::Vector{Tuple{Int,Int,Int,Float64}}
    dihedral::Vector{Tuple{Int,Int,Int,Int,Float64}}

    Embed2DConstraints() = new([],[],[],[])
end



function groupconstraint!(con::Embed2DConstraints, coords::Coordinate2DSubset)
    push!(con.groups, (keys, coords))
end
