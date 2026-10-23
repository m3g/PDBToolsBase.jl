"""
    center_of_mass(atoms::AbstractVector{<:Atom})

Calculate the center of mass of the atoms.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.test_file("small.pdb"));

julia> center_of_mass(atoms)
3-element StaticArraysCore.SVector{3, Float32} with indices SOneTo(3):
  -5.5844226
 -13.1104145
  -7.1399713
```

"""
function center_of_mass(atoms::AbstractVector{<:Atom})
    totmass = mass(atoms)
    center = SVector(
        sum(at.x * mass(at) for at in atoms),
        sum(at.y * mass(at) for at in atoms),
        sum(at.z * mass(at) for at in atoms)
    ) ./ totmass
    return center
end

@testitem "center_of_mass" setup = [PublicNames] begin
    using PDBToolsBase
    atoms = [Atom(name="C", x=-1.0, y=-1.0, z=-1.0), Atom(name="C", x=1.0, y=1.0, z=1.0)]
    @test center_of_mass(atoms) ≈ [0.0, 0.0, 0.0] atol = 1e-10
    atoms = read_pdb(PDBToolsBase.SMALLPDB)
    @test center_of_mass(atoms) ≈ [-5.584422752942997, -13.110413157869903, -7.139970815730879] atol = 1e-5
end

"""
    moveto!(atoms::AbstractVector{<:Atom}; center::AbstractVector{<:Real}=SVector(0.0, 0.0, 0.0))

Move the center of mass of the atoms to the specified `center` position, which defaults to the origin.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.test_file("small.pdb"));

julia> center_of_mass(atoms)
3-element StaticArraysCore.SVector{3, Float32} with indices SOneTo(3):
  -5.5844226
 -13.1104145
  -7.1399713

julia> moveto!(atoms; center = [1.0, 2.0, 3.0]);

julia> center_of_mass(atoms)
3-element StaticArraysCore.SVector{3, Float32} with indices SOneTo(3):
 1.0
 2.0000014
 3.0
```
"""
function moveto!(atoms::AbstractVector{<:Atom}; center::AbstractVector{<:Real}=SVector(0.0, 0.0, 0.0))
    cm = center_of_mass(atoms)
    for at in atoms
        at.x += center[1] - cm[1]
        at.y += center[2] - cm[2]
        at.z += center[3] - cm[3]
    end
    return atoms
end

@testitem "moveto!" setup = [PublicNames] begin
    using PDBToolsBase
    atoms = read_pdb(PDBToolsBase.SMALLPDB)
    @test center_of_mass(moveto!(atoms)) ≈ [0.0, 0.0, 0.0] atol = 1e-5
    @test center_of_mass(moveto!(atoms; center=[1.0, 2.0, 3.0])) ≈ [1.0, 2.0, 3.0] atol = 1e-5
end