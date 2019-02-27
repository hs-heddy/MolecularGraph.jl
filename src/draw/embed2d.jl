#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    compute2dcoords

using MolecularGraph.Geometry: _angle


"""
    compute2dcoords(mol::MolGraph) -> InternalCoordinates

Compute 2D coordinates of the molecular graph.
"""
function compute2dcoords(mol::VectorMol)
    # TODO: kekurize

    # Initial coordinates
    if is_outerplanar(mol)
        coords = outerplanar_embed2d(mol)
    else
        noprings = Set{Int}()
        nopnodes = Set{Int}()
        wd = sssrweakdual(mol)
        for bcon in biconnected_components(wd)
            for r in bcon
                push!(noprings, r)
                union!(nopnodes, mol.annotation[:Topology].rings[r])
            end
        end
        rset = Set(1:length(mol.annotation[:Topology].rings))
        oprings = setdiff(rset, noprings)
        opnodes = union(mol.annotation[:Topology].rings[r] for r in oprings)
        setdiff!(nopnodes, opnodes)
        # TODO: pending
        coords = outerplanar_embed2d(mol, pending=nopnodes)
        for r in noprings
            nodes = mol.annotation[:Topology].rings[r]
            # TODO: 2D embedding based on graph distance
            embedding = graphdistembedding(nodesubgraph(nodes))
            # TODO: assign graph dist emmbedding to the outerplanar embedding
            assign_cartesian2d!(coords, embedding)
        end
    end
    # so far, init coords and constraints (which nodes are fixed) are determined

    # TODO: apply chain geometry constraints (some rotations)
    # TODO: apply macrocycle constraints (add constraint and generate conformers)
    # TODO: generate conformers and evaluate
    # TODO: conformation optimization
    # TODO: force directed optimization
    return cartesian2d(coords)
end



struct OuterplanarEmbed2DState{G<:VectorMol}
    mol::G

    dfsindex::Dict{Int,Int} # node index, dfs tree index
    pred::Dict{Int,Int}
    geometry::Dict{Int,Symbol} # node index, :chain, :block, :spiro
    clockwise::Dict{Int,Bool} # node index, direction
    flip::Dict{Int,Bool} # node index, direction
    branchcount::Dict{Int,Int} # node index, number of branches
    branchindex::Dict{Int,Int} # node index, branch number

    coords::InternalCoords

    function OuterplanarEmbed2DState{G}(mol, root) where {G<:VectorMol}
        coords = internalcoords(nodecount(mol) + 3)
        # Set dummy nodes
        setcoord!(coords, 2, [1, nothing, nothing], [1.0, nothing, nothing])
        setcoord!(coords, 3, [2, 1, nothing], [1.0, 2 / 3, nothing])
        idx = Dict(-1 => 3, -2 => 2, -3 => 1)
        pred = Dict(-1 => -2, -2 => -3, -3 => -3, root => -1)
        geo = Dict(-1 => :chain, root => :chain)
        cw = Dict(-1 => true, root => false)
        flip = Dict(root => false)
        bcnt = Dict(-1 => 1)
        bidx = Dict(root => 1)
        new(mol, idx, pred, geo, cw, flip, bcnt, bidx, coords)
    end
end


"""
    outerplanar_embed2d(mol::VecterMol) -> Cartesian2D

Compute 2D embedding of the outerplanar graph which can be determined by a
simple DFS based algorithm.
"""
function outerplanar_embed2d(mol::G) where {G<:VectorMol}
    root = pop!(nodekeys(mol))
    state = OuterplanarEmbed2DState{G}(mol, root)
    dfs!(state, root)
    for (n, i) in state.dfsindex
        n > 0 && (state.coords.nodekeys[i] = n)
    end
    return state.coords
end


function dfs!(state::OuterplanarEmbed2DState, n::Int)
    state.dfsindex[n] = length(state.dfsindex) + 1
    p1 = state.pred[n]
    p2 = state.pred[p1]
    p3 = state.pred[p2]
    i_n = state.dfsindex[n]
    i1 = state.dfsindex[p1]
    i2 = state.dfsindex[p2]
    i3 = state.dfsindex[p3]
    i = get(state.branchindex, n, 0)
    bcnt = state.branchcount[p1]
    rmem = collect(get(state.mol[:RingMem], p1, Set{Int}()))
    ringsize = r -> length(state.mol.annotation[:Topology].rings[r])
    if length(rmem) == 0
        # chain -> chain
        if bcnt == 1
            # 1.0 is invalid (sp orbital cases)
            angle = get(state.mol[:Pi], p1, 0) == 2 ? 0.99 : 2 / 3
            dihedral = 1.0
        else
            a = 2.0i / (bcnt + 1)
            angle = a >= 1 ? a - 1.0 : a
            dihedral = a >= 1 ? 0.0 : 1.0
        end
    elseif length(rmem) == 1
        rsize = ringsize(rmem[1])
        ext = 1 / rsize
        if state.geometry[n] == :chain
            if state.geometry[p1] == :chain
                # block entry point -> chain
                bf = p2 > 0 ? i / (bcnt + 2) : (i - 1) / (bcnt + 1) # p2 < 0: root
                angle = bf
                dihedral = state.flip[n] ? 1.0 : 0.0
            else
                # block -> chain
                angle = ext + i / (bcnt + 1)
                dihedral = state.flip[n] ? 0.0 : 1.0
            end
        else
            if state.geometry[p1] == :chain
                # block entry point -> block
                bf = p2 > 0 ? 1 / (bcnt + 2) : 1 / (bcnt + 1) # p2 < 0: root
                angle = ext + bf
            else
                # block -> block
                angle = 1 - 2ext
            end
            dihedral = state.flip[n] ? 0.0 : 1.0
        end
    elseif length(rmem) == 2
        r1size = ringsize(rmem[1])
        r2size = ringsize(rmem[2])
        ext = 1 / r1size + 1 / r2size
        if state.geometry[n] == :chain
            if state.geometry[p1] == :chain
                # block entry point -> chain
                bf = p2 > 0 ? i / (bcnt + 2) : (i - 1) / (bcnt + 1) # p2 < 0: root
                angle = ext * 2 * bf
                dihedral = state.flip[n] ? 1.0 : 0.0
            else
                # block -> chain
                angle = ext * 2 * i / (bcnt + 1)
                dihedral = state.flip[n] ? 0.0 : 1.0
            end
        else
            # TODO: block entry point -> block
            # block -> block
            factor = Dict(:block => 2, :spiro => 1)
            angle = ext * factor[state.geometry[n]]
            dihedral = state.flip[n] ? 0.0 : 1.0
        end
    end
    setcoord!(state.coords, i_n, [i1, i2, i3], [1.0, angle, dihedral])

    # Node priority (block -> spiro -> chains)
    bin = Dict{Int,Set{Tuple{Int,Int}}}()
    block = Tuple{Int,Int}[]
    spiro = Tuple{Int,Int}[]
    chains = Tuple{Int,Int}[]
    for (nbr, bond) in neighbors(state.mol, n)
        if nbr in keys(state.dfsindex)
            continue # visited
        elseif length(state.mol[:RingBondMem][bond]) > 1
            continue # not Hamiltonian path
        elseif state.mol[:RingBond][bond]
            m = pop!(collect(state.mol[:RingBondMem][bond]))
            m in keys(bin) || (bin[m] = Set{Tuple{Int,Int}}())
            push!(bin[m], (nbr, bond))
        else
            push!(chains, (nbr, bond))
        end
    end
    for b in values(bin)
        if length(b) == 1
            push!(block, pop!(b))
        elseif length(b) == 2
            push!(spiro, pop!(b))
        end
    end

    # Block
    state.branchcount[n] = length(chains)
    if !isempty(block)
        (nbr, bond) = pop!(block)
        state.clockwise[nbr] = (state.geometry[n] == :chain
            || state.mol[:RingMemCount][n] == 2)
        state.flip[nbr] = state.clockwise[nbr] == state.clockwise[n]
        state.pred[nbr] = n
        state.geometry[nbr] = :block
        dfs!(state, nbr)
    end
    # Spiro
    if !isempty(spiro)
        (nbr, bond) = pop!(spiro)
        state.clockwise[nbr] = (state.geometry[n] == :chain
            || state.mol[:RingMemCount][n] == 2)
        state.flip[nbr] = state.clockwise[nbr] == state.clockwise[n]
        state.pred[nbr] = n
        state.geometry[nbr] = :spiro
        dfs!(state, nbr)
    end
    # Chains
    for (i, (nbr, bond)) in enumerate(chains)
        state.clockwise[nbr] = true
        state.flip[nbr] = (state.geometry[n] == :chain
            ? false : state.clockwise[nbr] == state.clockwise[n])
        state.pred[nbr] = n
        state.geometry[nbr] = :chain
        state.branchindex[nbr] = i
        dfs!(state, nbr)
    end

end


"""
    sssrweakdual(mol::VectorMol) -> MultiUDGraph

Compute a weak dual whose planar embedding represents SSSR.
"""
function sssrweakdual(mol::VectorMol)
    rcnt = length(state.mol.annotation[:Topology].rings)
    wd = MultiUDGraph(1:rcnt, [])
    ecnt = 0
    for mem in mol[:RingBondMem]
        length(mem) == 2 || continue
        memc = collect(mem)
        ecnt += 1
        updateedge!(wd, memc[1], memc[2], ecnt)
    end
    return wd
end


"""
    graphdistembedding(graph::UDGraph; kwargs...) -> Cartesian2D

Compute 2D embedding based on the graph distance.
"""
function graphdistembedding(graph::UDGraph)
    n = nodecount(mol)
    # TODO: adjacencymatrix
    D = adjacencymatrix(graph)
    H = LinearAlgebra.I - fill(1 / n, (n, n))
    G = - 0.5 * H * D * H # Gram matrix
    F = LinearAlgebra.eigen(G)
    od = sortperm(F.values)
    xidx = findall(i -> i == 1, od)
    yidx = findall(i -> i == 2, od)
    coords = zeros(n, 2)
    coords[:, 1] = F.vectors[:, xidx] * sqrt(F.values[xidx])
    coords[:, 2] = F.vectors[:, yidx] * sqrt(F.values[yidx])
    return Cartesian2D(coords)
end
