#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    forcedirected


mutable struct ForceConstraints
    mincycle::Dict{Int,Int} # node, minimum cycle size

    function ForceConstraints(cyclesize)
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
        state.maxiter = 50
        state.tick = 0.01
        state.stdlength = 1.0
        state.repulsion = 0.01
        state.spring = 2.0
        state.hinge = 0.002
        state.decay = 0.9
        return state
    end
end


"""
    forcedirected(graph::UDGraph) -> Cartesian2D

Compute 2D embedding based on the graph distance.
"""
function forcedirected(
        graph::G, init::Cartesian2D, constraints::ForceConstraints
        ) where {G<:UDGraph}
    state = ForceDirectedState{G}(graph, init, constraints)
    t = state.tick
    for i in 1:state.maxiter
        # Calculate forces
        f = zeros(size(state.coords))
        for n1 in nodekeys(graph)
            p1 = state.coords[n1, :]
            # Repulsion constraints
            for n2 in nodekeys(graph)
                n2 == n1 && continue
                p2 = state.coords[n2, :]
                factor = min(state.repulsion / norm(p2 - p1) ^ 2, state.repulsion * 100)
                f[n1, :] -= factor * normalize(p2 - p1)
            end
            # Spring constraint
            for nbr in neighborkeys(graph, n1)
                pnbr = state.coords[nbr, :]
                factor = (norm(pnbr - p1) - state.stdlength) * state.spring
                if degree(graph, n1) > 2
                    factor *= 2 # TODO This can work?
                end
                f[n1, :] += factor * normalize(pnbr - p1)
            end
            # Hinge constraint
            if degree(graph, n1) == 2
                (a, b) = neighborkeys(graph, n1)
                p1a = state.coords[a, :] - p1
                p1b = state.coords[b, :] - p1
                minc = state.constraints.mincycle[n]
                ang = cos(pi - 2pi / minc)
                cross = p1a[1] * p1b[2] - p1b[1] * p1a[2]
                direction = cross >= 0 ? -1 : 1
                factor = (dot(p1a, p1b) / (norm(p1a) * norm(p1b)) - ang) * state.hinge
                normveca = rotate(normalize(p1a - p1), pi / 2 * direction)
                f[a, :] += factor * normvec
                normvecb = rotate(normalize(p1b - p1), pi / 2 * direction)
                f[b, :] += factor * normvec
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
    end
    return Cartesian2D(state.coords)
end
