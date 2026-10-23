#
# Function names shared by the packages of the ecosystem (PDBTools, MolSimToolkit,
# ComplexMixtures, Packmol), which add methods to them. These were previously
# defined in MolSimToolkitShared.
#
function distance end
function positions end
function get_atoms end

"""
    distances

Shared function name: methods are defined by packages of the ecosystem (e.g. MolSimToolkit.jl),
and documented there.
"""
function distances end

"""
    bulk_coordination

Shared function name: methods are defined by packages of the ecosystem (e.g. MolSimToolkit.jl),
and documented there.
"""
function bulk_coordination end

"""
    save

Shared function name: methods are defined by packages of the ecosystem (e.g. PDBTools.jl,
ComplexMixtures.jl), and documented there.
"""
function save end

"""
    load

Shared function name: methods are defined by packages of the ecosystem (e.g. PDBTools.jl,
ComplexMixtures.jl), and documented there.
"""
function load end

@testitem "shared names" begin
    if VERSION >= v"1.11"
        @test Base.ispublic(PDBToolsBase, :distance)
        @test Base.ispublic(PDBToolsBase, :distances)
        @test Base.ispublic(PDBToolsBase, :bulk_coordination)
        @test Base.ispublic(PDBToolsBase, :positions)
        @test Base.ispublic(PDBToolsBase, :get_atoms)
        @test Base.ispublic(PDBToolsBase, :save)
        @test Base.ispublic(PDBToolsBase, :load)
    end
end

"""
    coordination_number

Shared function name: methods are defined by packages of the ecosystem (e.g. MolSimToolkit.jl,
ComplexMixtures.jl), and documented there. Calling it with arguments not supported by any
method throws an `ArgumentError`.
"""
function coordination_number(args...; kargs...)
    throw(ArgumentError("""\n
        Invalid arguments for the `coordination_number` function.
        Please check the documentation for the correct call signature, by typing:

        julia> ? coordination_number

    """))
end
@testitem "coordination_number" begin
    import PDBToolsBase: coordination_number
    # The coordination_number(args...; kargs...) is a placeholder for the docs only
    @test_throws ArgumentError coordination_number("string.txt", 1.0)
end
