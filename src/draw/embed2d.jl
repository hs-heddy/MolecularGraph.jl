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
function compute2dcoords(mol::VectorMol)
    # TODO: kekurize

    # Initial coordinates
    chains = Set{PointSet}()
    outerplanar = Set{PointSet}()
    nonouterplanar = Set{PointSet}()
    constraints = Embed2DConstraint()

    # Decompose
    for bcon in biconnected_components(mol)
        wd = sssrweakdual(bcon)
        for wbcon in biconnected_components(wd)
            nodes = Set{Int}()
            for r in bcon
                union!(nodes, mol.annotation[:Topology].rings[r])
            end
            push!(nonouterplanar, nodes)
        end
        push!(outerplanar, outerplanar_embed2d(mol, nodes))
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
