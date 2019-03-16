#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    decompose_chain,
    chain_embed


struct DecomposeChainState{G<:UndirectedGraph}
    graph::G
    biconn::Set{Int}
    chains::Set{Vector{Int}}

    queue::Vector{Int}
    pred::Dict{Int,Union{Int,Nothing}}

    function DecomposeChainState{G}(graph, root) where {G<:VectorMol}
        new(graph, biconn, Set())
    end
end


"""
    decompose_chain(graph::) -> Cartesian2D

Decompose chain subgraphs
"""
function decompose_chain(graph::G, biconn::Set{Int}
        ) where {G<:UndirectedGraph}
    chains = []
    for n in state.biconn
        state.pred[n] = nothing
        for nbr in neighborkeys(graph, n)
            nbr in state.biconn && continue # Bicoonnected component nodes
            nbr == state.pred[n] && continue # Already done
            dfs!(state, nbr)
            push!(chains, copy(state.queue))
            empty!(state.queue)
        end
    end
end


function dfs!(state::DecomposeChainState, n::Int)
    push!(state.queue, n)
    for (nbr, bond) in neighbors(state.graph, n)
        nbr == nothing && continue # Root
        nbr == state.pred[n] && continue # Predecessor
        nbr in state.biconn && continue # Bicoonnected component nodes
        dfs!(state, nbr)
    end
end


function chain_embed(state::DecomposeChainState)
    chains = Coordinate2DSubset[]
    for chain in state.chains
        nodes = Set(chain)
        while !isempty(nodes)
            subg = nodesubgraph(nodes)
            path = longestshortestpath(subg)
            coords = InternalCoordinates()
            flip = -1
            for p in path
                angle = state.constraint.angle[p] ? 2/3 : 179/180
                setcoords!(coords, n, [p1, p2, p3], [1, 2/3, flip])
                flip = flip == 1 ? -1 : 1
            end
            push!(chains, cartesian2dsubset(coords, path))
            setdiff!(nodes, Set(path))
        end
    end
end
