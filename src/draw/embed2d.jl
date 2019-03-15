#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    compute2dcoords,
    sssrweakdual


"""
    compute2dcoords(mol::MolGraph) -> InternalCoordinates

Compute 2D coordinates of the molecular graph.
"""
function compute2dcoords(mol::M) where {M<:VectorMol}
    # TODO: kekurize

    # Initial coordinates
    Embed2DConstraints()

    # Decompose
    bconstate = BiconnectedState{M}(mol)
    dfs!(bconstate)
    for conn in connected_components(edgesubgraph(bconstate.bridges))
        # Chains
        dec = decompose_chain(nodesubgraph(conn), bconstate)
        for chain in dec
            emb = chain_embed(nodesubgraph(chain), mol[:Pi])
            groupconstraint!(constraints, emb)
        end
    end
    for biconn in biconnected
        wd = sssrweakdual(nodesubgraph(biconn))
        wdbcs = biconnected_components(wd)
        # Non-outerplanar components
        for wdbcs in wdbcs
            nodes = Set{Int}()
            for r in wdbc
                union!(nodes, mol.annotation[:Topology].rings[r])
            end
            groupconstraint!(constraints, cartesian_embed(mol, nodes))
        end
        # Outerplanar components
        for r in setdiff(nodekeys(wd), union(wdbcs...))
            nodes = Set{Int}()
            union!(nodes, mol.annotation[:Topology].rings[r])
            groupconstraint!(constraints, outerplanar_embed(mol, nodes))
        end
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
