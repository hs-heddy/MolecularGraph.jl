#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    decompose


"""
    fragmentize(mol::MolGraph) -> Vector{Subgraph}

Decompose graph into chains and biconnected components
"""
function componentgraph(graph::UndirectedGraph)
    components = multigraph(1:circuitrank(graph), [])
    # Cycles
    ecnt = 0
    for (e, mem) in sssredgemem(graph)
        length(mem) < 2 && continue
        for (u, v) in combinations(mem)
            ecnt += 1
            edge = Edge(u, v, (e.u, e.v))
            updateedge!(components, edge, ecnt)
        end
    end
    # Adjacent to cycle
    cutvers = cutvertices(graph)
    while !isempty(cutvers)
        v = pop!(cutvers)
        if sssrmemcount(graph, v) == 2 # spiro
            ecnt += 1
            edge = Edge(sssrmem(graph, v)[1], sssrmem(graph, v)[2], (v,))
            updateedge!(components, edge, ecnt)
        elseif sssrmemcount(graph, v) == 1 # chain
            for nbr in neighborkeys(graph, v)
                state = DFSState(graph, v)
                dfs!(state, nbr)
                node = Node(state.path)
                setdiff!(cutvers, path)
                for t in state.term
                    edge = Edge(sssrmem(graph, v)[1], sssrmem(graph, v)[2], (v,))
                    updateedge!(components, edge, ecnt)
                end
                updatenode!(components, node, ncnt)
                updateedge!(components, edge, ecnt)
            end
        else
            continue
        end
    end
    # Isolated acyclic
    nocycle = setdiff(nodekeys(graph), biconnected_components(graph), )
    for comp in connected_components(nodesubgraph(nocycle))
        updatenode!(components, comp)
    end
end

function planarityminor!(graph::UndirectedGraph)
    for biconn in biconnected
        wd = sssrweakdual(nodesubgraph(biconn))
        # Non-outerplanar components
        nopset = Set{Int}()
        for wdbc in biconnected_components(wd)
            union!(nopset, wdbc)
            edges = union(cycleedges(mol, r) for r in wdbc)
            esub = edgesubgraph(mol, edges)
            push!(grouplinks, cartesian_embed(mol, nodekeys(esub)))
        end
        # Outerplanar components
        opset = setdiff(nodekeys(wd), nopset)
        opedges = union(cycleedges(mol, r) for r in opset)
        for conn in connected_components(edgesubgraph(opedges))
            push!(grouplinks, outerplanar_embed(mol, conn))
        end
        for edge in setdiff(edgekeys(wd), nopedges, opedges)
            links[edge.u][edge.v] = wd
            links[edge.v][edge.u] = wd
        end
    end


end


"""
    sssrweakdual(mol::VectorMol) -> MultiUDGraph

Compute a weak dual whose planar embedding represents SSSR.
"""
function sssrweakdual(graph::UndirectedGraph)
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


function edgecontraction(graph, nodes; allow_multi=true)
    newnodes = setdiff(nodekeys(graph), nodes)
    newg = newgraph(nodesubgraph(graph, newnodes))
    newi = nodecount(graph) + 1
    updatenodes!(newg, Node(), newi)
    for n in nodes
        for nbr in neighborkeys(graph, n)
            updateedges!(newg, Edge(nbr, newi))
        end
    end
    return newg
end
