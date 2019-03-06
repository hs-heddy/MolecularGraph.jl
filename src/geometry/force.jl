#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    forcedirected


mutable struct ForceConstraints
    mincycle::Dict{Int,Int} # node, minimum cycle size

    function ForceConstraints(cyclesize::Vector{Set{Int}})
        mincyc = Dict(i => maximum(c) for (i, c) in enumerate(cyclesize))
        new(mincyc)
    end
end


mutable struct ForceDirectedState{G<:UDGraph}
    graph::G
    constraints::ForceConstraints

    velo::Matrix{Float64}
    force::Matrix{Float64}
    maxiter::Int
    tick::Float64
    stdlength::Float64
    repulsion::Float64
    spring::Float64
    hinge::Float64
    decay::Float64

    coords::Matrix{Float64}

    function ForceDirectedState{G}(graph, init, constraints) where {G<:UDGraph}
        state = new()
        initmat = rawdata(init)
        state.coords = deepcopy(initmat)
        state.constraints = constraints
        state.velo = zeros(size(initmat))
        state.force = zeros(size(initmat))
        state.maxiter = 100
        state.tick = 0.01
        state.stdlength = 1.0
        state.repulsion = 1.0
        state.spring = 2.0
        state.hinge = 2.0
        state.decay = 1.0
        return state
    end
end


"""
    forcedirected(graph::UDGraph) -> Coordinates

Optimize cartesian embedding by force simulation.
"""
function forcedirected(
        graph::G, init::C, constraints::ForceConstraints
        ) where {G<:UDGraph, C<:Coordinates}
    state = ForceDirectedState{G}(graph, init, constraints)
    t = state.tick
    for i in 1:state.maxiter
        # Calculate forces
        f = zeros(size(state.coords))
        for n1 in nodekeys(graph)
            p1 = state.coords[n1, :]
            # Repulsion constraints
            for n2 in setdiff(nodekeys(graph), neighborkeys(graph, n1))
                n2 == n1 && continue
                p2 = state.coords[n2, :]
                factor = min(state.repulsion / norm(p2 - p1) ^ 2, state.repulsion * 10)
                f[n1, :] -= factor * normalize(p2 - p1)
            end
            # Spring constraints
            for nbr in neighborkeys(graph, n1)
                pnbr = state.coords[nbr, :]
                factor = state.spring * log(2, norm(pnbr - p1) / state.stdlength)
                f[n1, :] += factor * normalize(pnbr - p1)
            end
            # Hinge constraints
            if degree(graph, n1) == 2 && init isa Cartesian2D
                (a, b) = neighborkeys(graph, n1)
                p1a = state.coords[a, :] - p1
                p1b = state.coords[b, :] - p1
                minc = state.constraints.mincycle[n1]
                idealang = cos(pi - 2pi / minc)
                realang = dot(p1a, p1b) / (norm(p1a) * norm(p1b))
                torsion = realang - idealang
                cross = p1a[1] * p1b[2] - p1a[2] * p1b[1]
                adir = cross >= 0 ? -1 : 1
                bdir = cross >= 0 ? 1 : -1
                factor = state.hinge * torsion
                rotation = x -> [cos(x) -sin(x); sin(x) cos(x)]
                normveca = rotation(pi / 2 * adir) * normalize(p1a - p1)
                af = factor * normveca
                normvecb = rotation(pi / 2 * bdir) * normalize(p1b - p1)
                bf = factor * normvecb
                f[a, :] += af
                f[b, :] += bf
                f[n1, :] -= af + bf
            end
        end
        # Apply forces
        energy = 0.0
        for n1 in nodekeys(graph)
            state.velo[n1, :] += (state.force[n1, :] + f[n1, :]) / 2 * t
            state.force[n1, :] = f[n1, :]
            state.coords[n1, :] += state.velo[n1, :] * t + state.force[n1, :] * (t ^ 2) * 0.5
            energy += norm(state.force[n1, :])
        end
        println(energy)
        energy < 0.001 && break
        t *= state.decay
    end
    return C(state.coords)
end
