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
    componentgraph = componentgraph(mol)
    embed!(componentgraph)

    # TODO: apply macrocycle constraints (add constraint and generate conformers)
    # TODO: generate conformers and evaluate
    # TODO: conformation optimization
    # TODO: force directed optimization
    # TODO: place molecules
    return cartesian2d(coords)
end
