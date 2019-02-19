
@testset "draw.embed2d" begin

@testset "embed2d" begin
    disp!(mol, coords) = begin
        mol.coords[:Cartesian2D] = coords
        ASSETS_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets")
        dest = open(joinpath(ASSETS_DIR, "image", "test.svg"), "w")
        write(dest, drawsvg!(mol, 200, 200))
    end

    furan = smilestomol("NCC=1OC=CC1C")
    coords = compute2dcoords(furan)
    dist = Geometry.distance(segment(coords, 2, 8))
    @test isapprox(dist, 2.176, atol=1e-3)

    diazepam = smilestomol("CN1C2=C(C(C3=CC=CC=C3)=NCC1=O)C=C(Cl)C=C2")
    coords = compute2dcoords(diazepam)
    dist = Geometry.distance(segment(coords, 1, 15))
    @test isapprox(dist, 1.868, atol=1e-3)
    dist = Geometry.distance(segment(coords, 7, 16))
    @test isapprox(dist, 1.470, atol=1e-3)

    thiamine = smilestomol("OCCc1c(C)[n+](=cs1)Cc2cnc(C)nc(N)2")
    coords = compute2dcoords(thiamine)
    dist = Geometry.distance(segment(coords, 10, 15))
    @test isapprox(dist, 4.0, atol=1e-3)
    dist = Geometry.distance(segment(coords, 3, 10))
    @test isapprox(dist, 3.520, atol=1e-3)

    pyromellitimide = smilestomol("C1(=O)NC(=O)C=2C1=CC=3C(=O)NC(=O)C=3C=2")
    coords = compute2dcoords(pyromellitimide)
    dist = Geometry.distance(segment(coords, 2, 11))
    @test isapprox(dist, 4.252, atol=1e-3)
    dist = Geometry.distance(segment(coords, 5, 14))
    @test isapprox(dist, 4.252, atol=1e-3)

    spiro = smilestomol("C1CC12CC2")
    coords = compute2dcoords(spiro)
    dist = Geometry.distance(segment(coords, 1, 4))
    @test isapprox(dist, 1.732, atol=1e-3)
    dist = Geometry.distance(segment(coords, 1, 5))
    @test isapprox(dist, 2.0, atol=1e-3)
    # disp!(diazepam, coords)
    # display(coords)
end

end # draw.embed2d
