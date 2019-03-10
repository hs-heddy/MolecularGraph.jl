#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    compute2dcoords,
    outerplanar_embed2d,
    sssrweakdual


struct Embed2DConstraints
    groups::Vector{Tuple{Vector{Int}, Coordinates}}

    Embed2DConstraints() = new([])
end


function addgroupconstraint!(con::Embed2DConstraints, keys::Vector{Int},
                             coords::Coordinates)
    push!(con.groups, (keys, coords))
end


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
        nopcomps = Set{Set{Int}}()
        # Extract non-outerplanar components
        wd = sssrweakdual(mol)
        for bcon in biconnected_components(wd)
            nodes = Set{Int}()
            for r in bcon
                union!(nodes, mol.annotation[:Topology].rings[r])
            end
            push!(nopcomps, nodes)
        end
        # Graph distance based cartesian embedding
        econ = Embed2DConstraints()
        for nop in nopcomps
            subg = nodesubgraph(nop)
            graphem = graphdistembedding(subg)
            constraints = ForceConstraints(mol[:RingSize])
            if graphem isa Cartesian2D
                addgroupconstraint!(
                    econ, nodekeys(subg),
                    forcedirected2d(subg, graphem, constraints))
            elseif graphem isa Cartesian3D
                addgroupconstraint!(
                    econ, nodekeys(subg),
                    forcedirected3d(subg, graphem, constraints))
            end
        end
        # Combine outerplanar embedding and cartesian embedding
        coords = outerplanar_embed2d(mol, constraints=constraints)
    end
    # so far, init coords and constraints (which nodes are fixed) are determined

    # TODO: apply chain geometry constraints (some rotations)
    # TODO: apply macrocycle constraints (add constraint and generate conformers)
    # TODO: generate conformers and evaluate
    # TODO: conformation optimization
    # TODO: force directed optimization
    # TODO: place molecules
    return cartesian2d(coords)
end



struct OuterplanarEmbed2DState{G<:VectorMol}
    mol::G
    constraints::Embed2DConstraints

    dfsindex::Dict{Int,Int} # node index, dfs tree index
    pred::Dict{Int,Int}
    geometry::Dict{Int,Symbol} # node index, :chain, :block, :spiro
    clockwise::Dict{Int,Bool} # node index, direction
    flip::Dict{Int,Bool} # node index, direction
    branchcount::Dict{Int,Int} # node index, number of branches
    branchindex::Dict{Int,Int} # node index, branch number

    coords::InternalCoords

    function OuterplanarEmbed2DState{G}(mol, constraints, root) where {G<:VectorMol}
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
        new(mol, constraints, idx, pred, geo, cw, flip, bcnt, bidx, coords)
    end
end


"""
    outerplanar_embed2d(mol::VecterMol) -> Cartesian2D

Compute 2D embedding of the outerplanar graph which can be determined by a
simple DFS based algorithm.
"""
function outerplanar_embed2d(
        mol::G, constraints=Embed2DConstraints()::Embed2DConstraints
        ) where {G<:VectorMol}
    root = pop!(nodekeys(mol))
    state = OuterplanarEmbed2DState{G}(mol, constraints, root)
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
    # TODO: If it starts from ring node, block entry point should be more simple
    if state.geometry[n] == :fixed
        if state.geometry[p1] != :fixed
            angle =
        else
            angle = interiorangle(coords[p2-p1], coords[n-p1])
        end
        dihedral = state.flip[n] ? 1.0 : 0.0
    elseif state.geometry[p1] == :fixed
        
    elseif length(rmem) == 0
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

    # Node priority (block -> spiro -> fixed -> chains)
    bin = Dict{Int,Set{Int}}()
    pqueue = Tuple{Int,Symbol,Int}[]
    spiro = Tuple{Int,Symbol,Int}[]
    fixed = Tuple{Int,Symbol,Int}[]
    chains = Tuple{Int,Symbol,Int}[]
    chaincnt = 1
    for (nbr, bond) in neighbors(state.mol, n)
        if nbr in keys(state.dfsindex)
            continue # visited
        elseif length(state.mol[:RingBondMem][bond]) == 1 || (nbr in constraints.group
                && intersect(state.mol[:RingMem][n], state.mol[:RingMem][nbr]) > 1)
            m = pop!(state.mol[:RingBondMem][bond])
            m in keys(bin) || (bin[m] = Set{Tuple{Int,Int}}())
            push!(bin[m], nbr)
        elseif nbr in constraints.group
            push!(fixed, (nbr, :fixed, 0)) # Fixed coord
        elseif length(state.mol[:RingBondMem][bond]) > 1
            continue # not Hamiltonian path
        else
            push!(chains, (nbr, :chain, chaincnt))
            chaincnt += 1
        end
    end
    for b in values(bin)
        if length(b) == 1
            push!(pqueue, (pop!(b), :block, 0))
        elseif length(b) == 2
            push!(spiro, (pop!(b), :spiro, 0))
        end
    end
    append!(pqueue, spiro)
    append!(pqueue, chains)

    # Determine direction and move to next node
    for (nbr, geo, bidx) in pqueue
        state.pred[nbr] = n
        state.geometry[nbr] = :chain
        if geo == :chain
            state.clockwise[nbr] = true
            state.flip[nbr] = (state.geometry[n] == :chain
                ? false : state.clockwise[nbr] == state.clockwise[n])
            state.branchindex[nbr] = bidx
        elseif geo == :fixed
            state.clockwise[nbr] = (state.geometry[n] == :fixed
                ? cross(coords[n] - coords[p1], coords[nbr] - coords[n]) >= 0 : true)
            state.flip[nbr] = state.clockwise[nbr] == state.clockwise[n]
        else
            state.clockwise[nbr] = (state.geometry[n] == :chain
                || state.mol[:RingMemCount][n] == 2)
            state.flip[nbr] = state.clockwise[nbr] == state.clockwise[n]
        end
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
