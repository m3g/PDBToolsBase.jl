#
# Data for natural nucleosides
#
@kwdef struct NucleosideResidue
    name::String
    three_letter_code::String
    one_letter_code::String
    type::String
    mono_isotopic_mass::Float64
    mass::Float64
    charge::Int
    custom::Bool = false
end

#! format: off
const nucleoside_residues = OrderedDict{String,NucleosideResidue}(
    "ADO" => NucleosideResidue("Adenosine", "ADO", "A", "Purine",     267.094691, 267.094691, 0, false),
    "DA" => NucleosideResidue("Adenosine", "ADO", "A", "Purine",     267.094691, 267.094691, 0, false),
    "A" => NucleosideResidue("Adenosine", "ADO", "A", "Purine",     267.094691, 267.094691, 0, false),
    "GUO" => NucleosideResidue("Guanosine", "GUO", "G", "Purine",     283.090169, 283.090169, 0, false),
    "DG" => NucleosideResidue("Guanosine", "GUO", "G", "Purine",     283.090169, 283.090169, 0, false),
    "G" => NucleosideResidue("Guanosine", "GUO", "G", "Purine",     283.090169, 283.090169, 0, false),
    "URD" => NucleosideResidue("Uridine",   "URD", "U", "Pyrimidine", 244.068163, 244.068163, 0, false),
    "U" => NucleosideResidue("Uridine",   "URD", "U", "Pyrimidine", 244.068163, 244.068163, 0, false),
    "CYD" => NucleosideResidue("Cytidine",  "CYD", "C", "Pyrimidine", 244.068163, 244.068163, 0, false),
    "DC" => NucleosideResidue("Cytidine",  "CYD", "C", "Pyrimidine", 244.068163, 244.068163, 0, false),
    "C" => NucleosideResidue("Cytidine",  "CYD", "C", "Pyrimidine", 244.068163, 244.068163, 0, false),
    "THD" => NucleosideResidue("Thymidine", "THD", "T", "Pyrimidine", 244.068163, 244.068163, 0, false),
    "DT" => NucleosideResidue("Thymidine", "THD", "T", "Pyrimidine", 244.068163, 244.068163, 0, false),
    "T" => NucleosideResidue("Thymidine", "THD", "T", "Pyrimidine", 244.068163, 244.068163, 0, false),
    "INO" => NucleosideResidue("Inosine",   "INO", "I", "Purine",     268.090169, 268.090169, 0, false),
)
#! format: on

@doc """
    PDBTools.NucleosideResidue

Properties of a nucleoside: `name`, `three_letter_code`, `one_letter_code`, `type`
(`"Purine"` or `"Pyrimidine"`), `mono_isotopic_mass`, `mass`, `charge`, and `custom` 
(whether it was added with [`add_nucleoside_residue!`](@ref)).
""" NucleosideResidue

@doc """
    PDBTools.nucleoside_residues

Dictionary mapping residue names, as found in structure files, to the properties of the
nucleosides (of type `PDBTools.NucleosideResidue`). Custom nucleosides are added with
[`add_nucleoside_residue!`](@ref).
""" nucleoside_residues

"""
    add_nucleoside_residue!(resname::AbstractString, reference_residue::PDBTools.NucleosideResidue)

Add a custom nucleoside residue to the list of nucleoside residues. Returns the added
`NucleosideResidue`. Use `remove_custom_nucleoside_residues!()` to undo.

# Example

```jldoctest
julia> using PDBTools

julia> remove_custom_nucleoside_residues!();

julia> add_nucleoside_residue!("MYN", PDBTools.nucleoside_residues["ADO"])
PDBToolsBase.NucleosideResidue("MYN", "ADO", "A", "Purine", 267.094691, 267.094691, 0, true)

julia> isnucleoside(Atom(resname="MYN"))
true

julia> remove_custom_nucleoside_residues!(); # clean up
```

"""
function add_nucleoside_residue!(resname::AbstractString, reference_residue::PDBToolsBase.NucleosideResidue)
    if haskey(PDBToolsBase.nucleoside_residues, resname)
        @warn """\n
            Residue $resname already exists in the list of nucleoside residues. Overwriting.

        """ _file = nothing _line = nothing
    end
    PDBToolsBase.nucleoside_residues[resname] = PDBToolsBase.NucleosideResidue(
        resname,
        reference_residue.three_letter_code,
        reference_residue.one_letter_code,
        reference_residue.type,
        reference_residue.mono_isotopic_mass,
        reference_residue.mass,
        reference_residue.charge,
        true,
    )
    return PDBToolsBase.nucleoside_residues[resname]
end

"""
    remove_custom_nucleoside_residues!()

Remove all custom nucleoside residues from the nucleoside residue list.

# Example

```jldoctest
julia> using PDBTools

julia> remove_custom_nucleoside_residues!();

julia> add_nucleoside_residue!("MYN", PDBTools.nucleoside_residues["ADO"])
PDBToolsBase.NucleosideResidue("MYN", "ADO", "A", "Purine", 267.094691, 267.094691, 0, true)

julia> isnucleoside(Atom(resname="MYN"))
true

julia> remove_custom_nucleoside_residues!();

julia> isnucleoside(Atom(resname="MYN"))
false
```

"""
remove_custom_nucleoside_residues!() = filter!(r -> !last(r).custom, nucleoside_residues)

@testitem "nucleoside residues" setup = [PublicNames] begin
    using PDBToolsBase
    remove_custom_nucleoside_residues!()

    @test isnucleoside(Atom(resname="ADO"))
    @test isnucleoside(Atom(resname="GUO"))
    @test isnucleoside(Atom(resname="URD"))
    @test isnucleoside(Atom(resname="CYD"))
    @test isnucleoside(Atom(resname="THD"))
    @test isnucleoside(Atom(resname="INO"))
    @test !isnucleoside(Atom(resname="ALA"))

    @test ispurine(Atom(resname="ADO"))
    @test ispurine(Atom(resname="GUO"))
    @test ispurine(Atom(resname="INO"))
    @test !ispurine(Atom(resname="URD"))

    @test ispyrimidine(Atom(resname="URD"))
    @test ispyrimidine(Atom(resname="CYD"))
    @test ispyrimidine(Atom(resname="THD"))
    @test !ispyrimidine(Atom(resname="ADO"))

    # custom residue add/remove
    add_nucleoside_residue!("MYN", PDBToolsBase.nucleoside_residues["ADO"])
    @test isnucleoside(Atom(resname="MYN"))
    @test ispurine(Atom(resname="MYN"))
    remove_custom_nucleoside_residues!()
    @test !isnucleoside(Atom(resname="MYN"))

    # selections
    pdb = read_pdb(PDBToolsBase.test_file("structure.pdb"))
    @test length(filter(isnucleoside, pdb)) == 0  # no nucleosides in protein test structure

    s = read_pdb(PDBToolsBase.test_file("7YWH.cif"))
    @test length(s) == 6108
    @test length(filter(isnucleoside, s)) == 6108
    @test length(filter(ispyrimidine, s)) == 3400
    @test length(filter(ispurine, s)) == 2708
    @test residuename(resname(s[1])) == "Adenosine"
    r = collect(eachresidue(s))
    @test resname(r[1]) == "DA"
    @test residuename(r[1]) == "Adenosine"
end