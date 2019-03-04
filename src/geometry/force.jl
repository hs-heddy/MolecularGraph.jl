#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    forcedirected


mutable struct ForceDirectedState{G<:UDGraph}
    graph::G

    velo::Matrix{Float64}
    force::Matrix{Float64}
    maxiter::Int
    tick::Float64
    stdlength::Float64
    stdangle::Float64
    repulsion::Float64
    spring::Float64
    hinge::Float64
    decay::Float64

    coords::Matrix{Float64}

    function ForceDirectedState{G}(graph, init) where {G<:UDGraph}
        state = new()
        initmat = rawdata(init)
        state.coords = deepcopy(initmat)
        state.velo = zeros(size(initmat))
        state.force = zeros(size(initmat))
        state.maxiter = 50
        state.tick = 0.01
        state.stdlength = 1.0
        state.stdangle = -2/3
        state.repulsion = 0.01
        state.spring = 2.0
        state.hinge = 0.002
        state.decay = 0.9
        return state
    end
end


mutable struct ForceConstraints
    angle::Vector{Tuple{Int,Int,Int,Float64}} # center, end1, end2, cos(angle)

    function Constraints(graph, cycles, cycleedgesets)
        angles = Tuple{Int,Int,Int,Float64}[]
        for n in nodekeys(graph)
            for ((a, ea), (b, eb)) in combinations(neighbors(graph, n))
                c = intersect(cycleedgesets[ea], cycleedgesets[eb])
                size = length(cycles[pop!(c)])
                push!(consts, (n, a, b, cos(pi - 2pi / size))
            end
        end
        new(angles)
    end
end


"""
    forcedirected(graph::UDGraph) -> Cartesian2D

Compute 2D embedding based on the graph distance.
"""
function Geometry.forcedirected(
        graph::G, init::Cartesian2D, constraints::ForceConstraints
        ) where {G<:UDGraph}
    state = ForceDirectedState{G}(graph, init)
    t = state.tick
    for i in 1:state.maxiter
        energy = 0.0
        for n1 in nodekeys(graph)
            p1v = state.velo[n1, :] * t
            p1a = state.force[n1, :] * (t ^ 2) * 0.5
            state.coords[n1, :] +=  p1v + p1a
            p1 = state.coords[n1, :]
            force = [0.0, 0.0]
            # Repulsion constraints
            for n2 in nodekeys(graph)
                n2 == n1 && continue
                p2 = state.coords[n2, :]
                factor = min(state.repulsion / norm(p2 - p1) ^ 2, state.repulsion * 100)
                force -= factor * normalize(p2 - p1)
            end
            # Spring constraint
            for nbr in neighborkeys(graph, n1)
                pnbr = state.coords[nbr, :]
                factor = (norm(pnbr - p1) - state.stdlength) * state.spring
                force += factor * normalize(pnbr - p1)
            end
            # Hinge constraint
            for (a, b) in combinations(neighbors(graph, n1))
                pa = state.coords[a, :]
                pb = state.coords[b, :]
                factor = (interiorangle(pa, pb) - state.stdangle) * state.hinge
                force += factor * normalize(p2 - p1)
            end
            state.velo[n1, :] += (state.force[n1, :] + force) / 2 * t
            state.force[n1, :] = force
            energy += norm(state.force[n1, :])
        end
        println(energy)
        energy < 0.001 && break
    end
    return Cartesian2D(state.coords)
end
