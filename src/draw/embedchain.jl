#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    decompose_chain
    chain_embed


struct DecomposeChainState{G<:UndirectedGraph}
    graph::G
    biconn::Set{Int}
    chains::Set{Vector{Int}}

    function DecomposeChainState{G}(graph, root) where {G<:VectorMol}
        new(graph, biconn, Set())
    end
end


"""
    decompose_chain(graph::) -> Cartesian2D

Compute 2D embedding of the outerplanar graph which can be determined by a
simple DFS based algorithm.
"""
function decompose_chain(graph::G, biconn::Set{Int}
        ) where {G<:UndirectedGraph}

    state = DecomposeChainState{G}(graph, biconn)
    dfs!(state, pop!(nodekeys(graph)))

    chains = []
    for n in state.biconn
        nbrs = neighborkeys(graph, n)
    end
end


function dfs!(state::OuterplanarEmbedState, n::Int)
    for (nbr, bond) in neighbors(state.graph, n)
        nbr == state.pred[n] && continue # Predecessor
        state.ringbondcount[bond] == 2 && continue # Not Hamiltonian path
        state.pred[nbr] = n
        state.clockwise[nbr] = state.ringcount[n] == 2
        state.flip[nbr] = state.clockwise[nbr] == state.clockwise[n]
        dfs!(state, nbr)
    end
end
