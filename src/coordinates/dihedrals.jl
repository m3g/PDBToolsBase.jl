
"""
    dihedral(at1::Atom, at2::Atom, at3::Atom, at4::Atom)

Computes the dihedral angle given four atoms of type PDBTools.Atom.

# Example

```jldoctest; filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> using PDBTools

julia> pdb = read_pdb(PDBTools.test_file("structure.pdb"));

julia> C1 = pdb[11]; N2 = pdb[13]; CA2 = pdb[15]; C2 = pdb[22];

julia> phi = dihedral(C1, N2, CA2, C2) 
-36.70359f0
```

"""
dihedral(at1::Atom, at2::Atom, at3::Atom, at4::Atom) =
    dihedral(position(at1), position(at2), position(at3), position(at4))

@testitem "dihedral" setup = [PublicNames] begin
    using PDBToolsBase
    pdb = read_pdb(PDBToolsBase.test_file("structure.pdb"));
    N1 = pdb[1]; CA1 = pdb[5]; C1 = pdb[11]; N2 = pdb[13];
    @test dihedral(N1, CA1, C1, N2) ≈ 64.07296f0
    C1 = pdb[11]; N2 = pdb[13]; CA2 = pdb[15]; C2 = pdb[22];
    @test dihedral(C1, N2, CA2, C2) ≈ -36.70359f0
end
