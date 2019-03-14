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



function addgroupconstraint!(con::Embed2DConstraints, keys::Vector{Int},
                             coords::Coordinates)
    push!(con.groups, (keys, coords))
end
