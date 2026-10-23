"""
    Residue

Residue data structure. 

The Residue structure carries the properties of the residue or molecule of the atoms
it contains, but it does not copy the original vector of atoms, only the residue
meta data for each residue. Thus, changes in the residue atoms will be reflected in the
original vector of atoms.

### Example

```jldoctest
julia> using PDBTools

julia> pdb = wget("1LBD");

julia> residues = collect(eachresidue(pdb))
238-element Vector{Residue}[
    SER225A
    ALA226A
    ⋮
    MET461A
    THR462A
]

julia> resnum.(residues[1:3])
3-element Vector{Int32}:
 225
 226
 227

julia> residues[5].chain
"A"

julia> residues[8].range
52:58

julia> mass(residues[1])
82.0385f0

```

"""
struct Residue{T<:Atom,Vec<:AbstractVector{T}} <: AbstractStructuralElement{T}
    atoms::Vec
    range::UnitRange{Int}
    resname::StringType
    residue::Int32
    resnum::Int32
    chain::StringType
    model::Int32
    segname::StringType
end

# Necessary for the interface: define the same_struct_element function
same_struct_element(::Type{Residue}, at1::Atom, at2::Atom) = same_residue(at1, at2)
Base.in(at::Atom, r::Residue) = same_struct_element(Residue, at, first(r))

# Constructors
function Residue(atoms::AbstractVector{<:Atom}, range::AbstractRange{<:Integer})
    _check_unique(Residue, atoms, range)
    at1 = first(atoms[range])
    Residue(atoms, UnitRange{Int}(range), resname(at1), residue(at1), resnum(at1), chain(at1), model(at1), segname(at1))
end
Residue(atoms::AbstractVector{<:Atom}) = Residue(atoms, 1:length(atoms))


"""
    eachresidue(atoms::AbstractVector{<:Atom})

Iterator for the residues (or molecules) of a selection. 

### Example

```jldoctest
julia> using PDBTools

julia> atoms = wget("1LBD");

julia> eachresidue(atoms)
 Residue iterator with length = 238

julia> collect(eachresidue(atoms))
238-element Vector{Residue}[
    SER225A
    ALA226A
    ⋮
    MET461A
    THR462A
]
```

"""
eachresidue(atoms::AbstractVector{<:Atom}) = EachStructuralElement{Residue}(atoms)

# Specific getters for this type

"""
    name(residue::Residue)

Returns the name (residue name) of the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `String`: The residue name.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> name(residues[1])
"ALA"
```

"""
name(residue::Residue) = residue.resname

"""
    resname(residue::Residue)

Returns the residue name of the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `String`: The residue name.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> resname(residues[1])
"ALA"
```

"""
resname(residue::Residue) = residue.resname

"""
    residuename(r::Residue)

Returns the full residue name of the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `String`: The full residue name.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> residuename(residues[1])
"Alanine"
```

"""
residuename(r::Residue) = residuename(resname(r))

"""
    residue(residue::Residue)

Returns the sequential residue (molecule) number of the residue in the structure.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `Int32`: The sequential residue number.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> residue(residues[1])
1
```

"""
residue(residue::Residue) = residue.residue

"""
    resnum(residue::Residue)

Returns the residue number as written in the PDB file.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `Int32`: The residue number from the PDB file.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> resnum(residues[1])
1
```

"""
resnum(residue::Residue) = residue.resnum

"""
    chain(residue::Residue)

Returns the chain identifier of the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `String`: The chain identifier.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> chain(residues[1])
"A"
```

"""
chain(residue::Residue) = residue.chain

"""
    model(residue::Residue)

Returns the model number of the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `Int32`: The model number.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> model(residues[1])
1
```

"""
model(residue::Residue) = residue.model

"""
    segname(residue::Residue)

Returns the segment name of the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `String`: The segment name.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> segname(residues[1])
"PROT"
```

"""
segname(residue::Residue) = residue.segname

"""
    mass(residue::Residue)

Returns the total mass of all atoms in the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `Float32`: The total mass in atomic mass units (amu).

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> mass(residues[1])
73.09489f0
```

"""
mass(residue::Residue) = mass(@view residue.atoms[residue.range])

"""
    get_atoms(residue::Residue)

Returns a view of the atoms that belong to the residue.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `SubArray{Atom}`: A view of the atoms in the residue.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.SMALLPDB);

julia> residues = collect(eachresidue(atoms));

julia> res_atoms = get_atoms(residues[1]);

julia> length(res_atoms)
12
```

"""
get_atoms(residue::Residue) = @view(residue.atoms[residue.range])

"""
    charge(residue::Residue)

Returns the formal charge of the residue. For protein residues, the charge is obtained
from the standard protein residue database. For non-protein residues, a warning is issued
and 0 is returned.

# Positional Arguments

- `residue::Residue`: A `Residue` object.

# Returns

- `Int`: The formal charge of the residue.

# Example

```jldoctest
julia> using PDBTools

julia> atoms = read_pdb(PDBTools.TESTPDB, "protein and resname GLU");

julia> residues = collect(eachresidue(atoms));

julia> charge(residues[1])
-1
```

"""
function charge(residue::Residue)  
    if residue.resname in keys(protein_residues)
        return protein_residues[residue.resname].charge
    else
        @warn "Residue $(residue.resname) not found in protein residues. Assuming residue charge = 0"
        return 0
    end
end

@testitem "Residue iterator" setup = [PublicNames] begin
    using PDBToolsBase
    atoms = read_pdb(PDBToolsBase.test_file("structure.pdb"), isprotein)
    residues = eachresidue(atoms)
    @test length(residues) == 104
    @test name(Residue(atoms, 1:12)) == "ALA"
    @test resname(Residue(atoms, 1:12)) == "ALA"
    @test Residue(atoms, 1:12).range == 1:12
    @test Residue(atoms[1:12]).range == 1:12
    @test firstindex(residues) == 1
    @test lastindex(residues) == 104
    @test resname(last(residues)) == "THR"
    @test residue(last(residues)) == 104
    @test chain(last(residues)) == "A"
    @test_throws ArgumentError residues[1]
    residues = collect(eachresidue(atoms))
    @test atoms[1] in residues[1]
    @test !(atoms[1] in residues[2])
    @test atoms[1458] in residues[104]
    @test index.(filter(at -> name(at) in ("N", "HG1"), residues[2])) == [13, 21]
    @test length(get_atoms(residues[1])) == 12
    @test findall(at -> name(at) in ("N", "HG1"), residues[2]) == [1, 9]
    @test_throws ArgumentError Residue(atoms, 1:15)
end

#
# io show functions
#
function Base.show(io::IO, residue::Residue)
    compact = get(io, :compact, false)::Bool
    if compact
        print(io, "$(name(residue))$(resnum(residue))$(chain(residue))")
    else
        println(io, " Residue of name $(name(residue)) with $(length(residue)) atoms.")
        show(IOContext(io, :type => false), @view residue.atoms[residue.range])
    end
end

@testitem "residue show" setup = [PublicNames] begin
    using PDBToolsBase
    using ShowMethodTesting
    ENV["LINES"] = 10
    ENV["COLUMNS"] = 120
    ats = read_pdb(PDBToolsBase.SMALLPDB)
    r = eachresidue(ats)
    @test parse_show(r; repl=Dict("PDBToolsBase." => "")) ≈ "Residue iterator with length = 3"
    rc = collect(r)
    @test parse_show(rc; repl=Dict("PDBToolsBase." => "")) ≈ """
    3-element Vector{Residue}[ 
        ALA1A
        CYS2A
        ASP3A
    ]
    """
    @test parse_show(rc[1]; repl=Dict(r"^((?:[^\n]*\n){3}).*"s => s"\1")) ≈ """
     Residue of name ALA with 12 atoms.
    index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    PROT         1
    """
end

#
# Properties of residues
#
# Deal with alternate conformers with the second test
#
isprotein(residue::Residue) =
    haskey(protein_residues, resname(residue)) ||
    (length(resname(residue)) == 4 && haskey(protein_residues, @view(resname(residue)[2:4])))

isnucleoside(residue::Residue) = haskey(nucleoside_residues, resname(residue))

isacidic(r::Residue) = isprotein(r) && protein_residues[r.resname].type == "Acidic"
isaliphatic(r::Residue) = isprotein(r) && protein_residues[r.resname].type == "Aliphatic"
isaromatic(r::Residue) = isprotein(r) && protein_residues[r.resname].type == "Aromatic"
isbasic(r::Residue) = isprotein(r) && protein_residues[r.resname].type == "Basic"
ischarged(r::Residue) = isprotein(r) && protein_residues[r.resname].charge != 0
isneutral(r::Residue) = isprotein(r) && protein_residues[r.resname].charge == 0
ishydrophobic(r::Residue) = isprotein(r) && protein_residues[r.resname].hydrophobic
ispolar(r::Residue) = isprotein(r) && protein_residues[r.resname].polar
isnonpolar(r::Residue) = isprotein(r) && !ispolar(r)

isacidic(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].type == "Acidic"
isaliphatic(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].type == "Aliphatic"
isaromatic(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].type == "Aromatic"
isbasic(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].type == "Basic"
ischarged(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].charge != 0
isneutral(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].charge == 0
ishydrophobic(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].hydrophobic
ispolar(atom::Atom) = isprotein(atom) && protein_residues[atom.resname].polar
isnonpolar(atom::Atom) = isprotein(atom) && !ispolar(atom)

ispurine(r::Residue) = isnucleoside(r) && nucleoside_residues[r.resname].type == "Purine"
ispyrimidine(r::Residue) = isnucleoside(r) && nucleoside_residues[r.resname].type == "Pyrimidine"
ispurine(atom::Atom) = isnucleoside(atom) && nucleoside_residues[atom.resname].type == "Purine"
ispyrimidine(atom::Atom) = isnucleoside(atom) && nucleoside_residues[atom.resname].type == "Pyrimidine"

const water_residues = ["HOH", "OH2", "TIP3", "TIP3P", "TIP4P", "TIP5P", "TIP7P", "SPC", "SPCE"]
iswater(r::Residue; water_residues=water_residues) = r.resname in water_residues
iswater(atom::Atom; water_residues=water_residues) = atom.resname in water_residues

#
# Documentation of the classification functions of atoms and residues
#
@doc """
    isprotein(atom::Atom)
    isprotein(residue::Residue)

Returns `true` if the atom belongs to, or the residue is, a protein residue (as defined in
`PDBTools.protein_residues`, see [`add_protein_residue!`](@ref)).
""" isprotein

@doc """
    isnucleoside(atom::Atom)
    isnucleoside(residue::Residue)

Returns `true` if the atom belongs to, or the residue is, a nucleoside (as defined in
`PDBTools.nucleoside_residues`, see [`add_nucleoside_residue!`](@ref)).
""" isnucleoside

@doc """
    isbackbone(atom::Atom; backbone_atoms=["N", "CA", "C", "O"])

Returns `true` if the atom is a protein backbone atom, that is, if it is a protein atom with
one of the names in `backbone_atoms`.
""" isbackbone

@doc """
    issidechain(atom::Atom; not_side_chain_atoms=["N", "CA", "C", "O", "HN", "H", "HA", "HT1", "HT2", "HT3"])

Returns `true` if the atom is a protein side chain atom, that is, if it is a protein atom 
with a name that is not in `not_side_chain_atoms`.
""" issidechain

@doc """
    iswater(atom::Atom; water_residues=["HOH", "OH2", "TIP3", "TIP3P", "TIP4P", "TIP5P", "TIP7P", "SPC", "SPCE"])
    iswater(residue::Residue; water_residues=...)

Returns `true` if the atom belongs to, or the residue is, a water molecule, that is, if the
residue name is one of `water_residues`.
""" iswater

@doc """
    isacidic(atom::Atom)
    isacidic(residue::Residue)

Returns `true` if the atom belongs to, or the residue is, an acidic protein residue.
""" isacidic

@doc """
    isaliphatic(atom::Atom)
    isaliphatic(residue::Residue)

Returns `true` if the atom belongs to, or the residue is, an aliphatic protein residue.
""" isaliphatic

@doc """
    isaromatic(atom::Atom)
    isaromatic(residue::Residue)

Returns `true` if the atom belongs to, or the residue is, an aromatic protein residue.
""" isaromatic

@doc """
    isbasic(atom::Atom)
    isbasic(residue::Residue)

Returns `true` if the atom belongs to, or the residue is, a basic protein residue.
""" isbasic

@doc """
    ischarged(atom::Atom)
    ischarged(residue::Residue)

Returns `true` if the atom belongs to, or the residue is, a charged protein residue.
""" ischarged

@doc """
    isneutral(atom::Atom)
    isneutral(residue::Residue)

Returns `true` if the atom belongs to, or the residue is, a neutral protein residue.
""" isneutral

@doc """
    ishydrophobic(atom::Atom)
    ishydrophobic(residue::Residue)

Returns `true` if the atom belongs to, or the residue is, a hydrophobic protein residue.
""" ishydrophobic

@doc """
    ispolar(atom::Atom)
    ispolar(residue::Residue)

Returns `true` if the atom belongs to, or the residue is, a polar protein residue.
""" ispolar

@doc """
    isnonpolar(atom::Atom)
    isnonpolar(residue::Residue)

Returns `true` if the atom belongs to, or the residue is, a non-polar protein residue.
""" isnonpolar

@doc """
    ispurine(atom::Atom)
    ispurine(residue::Residue)

Returns `true` if the atom belongs to, or the residue is, a purine nucleoside.
""" ispurine

@doc """
    ispyrimidine(atom::Atom)
    ispyrimidine(residue::Residue)

Returns `true` if the atom belongs to, or the residue is, a pyrimidine nucleoside.
""" ispyrimidine

@testitem "residue of atom" setup = [PublicNames] begin
    pdb = read_pdb(PDBToolsBase.test_file("structure.pdb"))
    glu = filter(at -> resname(at) == "GLU", pdb)
    @test isacidic(glu[1])
    @test !isaliphatic(glu[1])
    @test !isaromatic(glu[1])
    @test !isbasic(glu[1])
    @test ischarged(glu[1])
    @test !isneutral(glu[1])
    @test !ishydrophobic(glu[1])
    @test ispolar(glu[1])
    @test !isnonpolar(glu[1])
    @test !iswater(glu[1])
    phe = filter(at -> resname(at) == "PHE", pdb)
    @test !isacidic(phe[1])
    @test !isaliphatic(phe[1])
    @test isaromatic(phe[1])
    @test !isbasic(phe[1])
    @test !ischarged(phe[1])
    @test isneutral(phe[1])
    @test ishydrophobic(phe[1])
    @test !ispolar(phe[1])
    @test isnonpolar(phe[1])
    @test !iswater(glu[1])
    wat = filter(iswater, pdb)
    @test iswater(wat[1])
    pdb = [ Atom(resname="AGLY"), Atom(resname="BGLY")]
    @test length(eachresidue(pdb)) == 1
    @test resname.(eachresidue(pdb)) == ["AGLY"] 
end

@testitem "full residue" setup = [PublicNames] begin
    pdb = read_pdb(PDBToolsBase.test_file("structure.pdb"))
    glu_atoms = filter(at -> resname(at) == "GLU", pdb)
    glu = collect(eachresidue(glu_atoms))
    @test charge(glu[1]) == -1
    @test isacidic(glu[1])
    @test !isaliphatic(glu[1])
    @test !isaromatic(glu[1])
    @test !isbasic(glu[1])
    @test ischarged(glu[1])
    @test !isneutral(glu[1])
    @test !ishydrophobic(glu[1])
    @test ispolar(glu[1])
    @test !isnonpolar(glu[1])
    @test !iswater(glu[1])
    @test mass(glu[1]) ≈ 128.1077 atol = 1e-3
    @test mass(glu[begin]) ≈ mass(glu[end])
    @test segname(glu[begin]) == "PROT"
    @test segname(glu[end]) == "PROT"
    phe_atoms = filter(at -> resname(at) == "PHE", pdb)
    phe = collect(eachresidue(phe_atoms))
    @test charge(phe[1]) == 0
    @test !isacidic(phe[1])
    @test !isaliphatic(phe[1])
    @test isaromatic(phe[1])
    @test !isbasic(phe[1])
    @test !ischarged(phe[1])
    @test isneutral(phe[1])
    @test ishydrophobic(phe[1])
    @test !ispolar(phe[1])
    @test isnonpolar(phe[1])
    @test !iswater(glu[1])
    water = filter(at -> iswater(at) && resnum(at) <= 3, pdb)
    watresiter = eachresidue(water)
    watres = collect(eachresidue(water))
    @test iswater(watres[1])
    @test charge(watres[1]) == 0
end
