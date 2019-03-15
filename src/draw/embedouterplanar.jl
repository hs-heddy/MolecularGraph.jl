#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    outerplanar_embed


struct OuterplanarEmbedState{G<:VectorMol}
    graph::G

    ringcount::Dict{Int,Int} # node index, number of cycles
    ringbondcount::Dict{Int,Int} # edge index, number of cycles
    ringsize::Dict{Int,Vector{Int}} # edge index, list of cycle sizes
    path::Vector{Int} # Hamiltonian path
    dfsindex::Dict{Int,Int} # node index, dfs tree index
    pred::Dict{Int,Int}
    clockwise::Dict{Int,Bool} # node index, direction
    flip::Dict{Int,Bool} # node index, direction

    coords::InternalCoords

    function OuterplanarEmbedState{G}(graph, root) where {G<:VectorMol}
        coords = internalcoords(nodecount(mol) + 3)
        # Set dummy nodes
        setcoord!(coords, 2, [1, nothing, nothing], [1.0, nothing, nothing])
        setcoord!(coords, 3, [2, 1, nothing], [1.0, 2 / 3, nothing])
        idx = Dict(-1 => 3, -2 => 2, -3 => 1)
        pred = Dict(-1 => -2, -2 => -3, -3 => -3, root => -1)
        cw = Dict(-1 => true, root => false)
        flip = Dict(root => false)
        new(graph, rcnt, rbond, rsize, idx, pred, cw, flip, coords)
    end
end


"""
    outerplanar_embed(mol::VecterMol) -> Cartesian2D

Compute 2D embedding of the outerplanar graph which can be determined by a
simple DFS based algorithm.
"""
function outerplanar_embed(graph::G) where {G<:UndirectedGraph}
    root = pop!(nodekeys(graph))
    state = OuterplanarEmbedState{G}(graph, root)
    dfs!(state, root)
    for n in state.path
        coords!(state, n)
    end
    for (n, i) in state.dfsindex
        n > 0 && (state.coords.nodekeys[i] = n)
    end
    return cartesian2d(state.coords)
end


function dfs!(state::OuterplanarEmbedState, n::Int)
    push!(state.path, n)
    state.dfsindex[n] = length(state.path)
    for (nbr, bond) in neighbors(state.graph, n)
        nbr == state.pred[n] && continue # Predecessor
        state.ringbondcount[bond] == 2 && continue # Not Hamiltonian path
        state.pred[nbr] = n
        state.clockwise[nbr] = state.ringcount[n] == 2
        state.flip[nbr] = state.clockwise[nbr] == state.clockwise[n]
        dfs!(state, nbr)
    end
end


function coords!(state::OuterplanarEmbedState, n::Int)
    if length(state.ringsize[n]) == 1
        rsize = state.ringsize[n][1]
        angle = 1 - 2 / rsize
        dihedral = state.flip[n] ? 0.0 : 1.0
    else
        r1size = state.ringsize[n][1]
        r2size = state.ringsize[n][2]
        angle = 2 / r1size + 2 / r2size
        dihedral = state.flip[n] ? 0.0 : 1.0
    end
    p1 = state.pred[n]
    p2 = state.pred[p1]
    p3 = state.pred[p2]
    i_n = state.dfsindex[n]
    i1 = state.dfsindex[p1]
    i2 = state.dfsindex[p2]
    i3 = state.dfsindex[p3]
    setcoord!(state.coords, i_n, [i1, i2, i3], [1.0, angle, dihedral])
end
