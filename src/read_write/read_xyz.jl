"""
    read_xyz(xyz_file::AbstractString, selection_function::Function = all)
    read_xyz(xyz_data::IOBuffer, selection_function::Function = all)

Reads a file in the xyz format and returns a `Structure`.

The xyz format consists of one or more frames, each with:

- a line with the number of atoms,
- a comment line, which is stored as the `title` of the structure (from the first frame),
- one line per atom, with the element symbol (or atom name, or atomic number) and the
  `x`, `y`, and `z` coordinates. Additional columns are ignored.

The first column is stored, as is, in the `name` field of the atoms (and is written back
by [`write_xyz`](@ref)). The element symbol, if recognized from the first column, is stored
in the `pdb_element` field (atomic numbers are converted to element symbols).
Each frame is read as a different model (`model` field), and the atoms of each frame
form a single residue.

If a selection function is provided, only the atoms for which `selection_function(atom)`
is true will be read.

### Example

```julia-repl
julia> using PDBTools

julia> s = read_xyz(PDBTools.test_file("ethanol.xyz"))
   Structure{Atom{Nothing}} with 9 atoms
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    C     XXX     X        1        1   -0.748   -0.015    0.024  0.00  0.00     1       X         1
       2    C     XXX     X        1        1    0.558    0.420   -0.619  0.00  0.00     1       X         2
⋮
       8    H     XXX     X        1        1    0.570    1.511   -0.745  0.00  0.00     1       X         8
       9    H     XXX     X        1        1    1.510    0.138   -2.348  0.00  0.00     1       X         9

julia> s.title
"Ethanol"
```

"""
function read_xyz end

read_xyz(xyzdata::IOBuffer, selection_function::F=all) where {F<:Function} =
    _parse_xyz(xyzdata, selection_function)

function read_xyz(filename::AbstractString, selection_function::F=all) where {F<:Function}
    atoms = open(expanduser(filename), "r") do f
        _parse_xyz(f, selection_function)
    end
    return atoms
end

function _parse_xyz(xyzdata::Union{IOStream,IOBuffer}, selection_function::F) where {F<:Function}
    atoms = Atom{Nothing}[]
    lastatom = Atom{Nothing}()
    title = ""
    imodel = 0
    # Buffer for the fields of each line, reused for all atoms
    fields = SubString{String}[]
    while !eof(xyzdata)
        line = readline(xyzdata)
        isempty(strip(line)) && continue # blank lines between or after frames
        natoms = tryparse(Int, strip(line))
        if isnothing(natoms) || natoms < 0
            throw(ArgumentError("Invalid xyz file: expected the number of atoms, got: \"$line\""))
        end
        imodel += 1
        comment = readline(xyzdata)
        imodel == 1 && (title = String(strip(comment)))
        for iatom in 1:natoms
            if eof(xyzdata)
                throw(ArgumentError("Invalid xyz file: expected $natoms atoms in frame $imodel, found $(iatom - 1)."))
            end
            atom = read_atom_xyz(readline(xyzdata), fields, lastatom, iatom, imodel)
            selection_function(atom) && push!(atoms, atom)
            lastatom = atom
        end
    end
    seekstart(xyzdata)
    if length(atoms) == 0
        throw(ArgumentError("""\n
            Could not find any atom in xyz file matching the selection.

        """))
    end
    return Structure(atoms, MetaData(nothing, nothing, isempty(title) ? nothing : title))
end

function read_atom_xyz(
    record::String,
    fields::Vector{SubString{String}},
    lastatom::Atom,
    iatom::Int,
    imodel::Int,
)
    empty!(fields)
    for field in eachsplit(record)
        push!(fields, field)
    end
    if length(fields) < 4
        throw(ArgumentError("Invalid xyz file: could not read atom from line: \"$record\""))
    end
    atom = Atom{Nothing}(;
        index=index(lastatom) + 1,
        index_pdb=iatom,
        resnum=1,
        residue=imodel,
        model=imodel,
        x=_parse(Float32, fields[2]),
        y=_parse(Float32, fields[3]),
        z=_parse(Float32, fields[4]),
    )
    atom.name = _parse(StringType, fields[1]; alt=_alt(StringType))
    # The first column can be the atomic number
    atomic_number = tryparse(Int, fields[1])
    atom.pdb_element = if isnothing(atomic_number)
        something(element(atom), StringType("X"))
    else
        _element_symbol_from_atomic_number(atomic_number)
    end
    return atom
end

function _element_symbol_from_atomic_number(atomic_number::Int)
    for el in values(elements)
        el.atomic_number == atomic_number && !el.custom && return el.symbol_string
    end
    throw(ArgumentError("Invalid xyz file: unknown atomic number $atomic_number."))
end

@testitem "read_xyz" setup = [PublicNames] begin
    using PDBToolsBase
    s = read_xyz(PDBToolsBase.test_file("ethanol.xyz"))
    @test s isa Structure{Atom{Nothing}}
    @test length(s) == 9
    @test s.title == "Ethanol"
    @test isnothing(s.unitcell)
    @test isnothing(s.connectivity)
    @test name.(s) == ["C", "C", "O", "H", "H", "H", "H", "H", "H"]
    @test element.(s) == name.(s)
    @test index.(s) == 1:9
    @test index_pdb.(s) == 1:9
    @test all(==(1), resnum.(s)) && all(==(1), residue.(s)) && all(==(1), model.(s))
    @test (s[3].x, s[3].y, s[3].z) == (0.716f0, -0.191f0, -1.871f0)
    @test length(eachresidue(s)) == 1
    @test read_xyz(IOBuffer(read(PDBToolsBase.test_file("ethanol.xyz"), String))) == s
    # coordinates are the same as those of the PDB file
    pdb = read_pdb(PDBToolsBase.test_file("ethanol.pdb"))
    @test positions(s) == positions(pdb)
    @test element.(s) == element.(pdb)
    # selection
    s = read_xyz(PDBToolsBase.test_file("ethanol.xyz"), at -> element(at) != "H")
    @test length(s) == 3
    @test index.(s) == 1:3
    # multiple frames, atom labels, atomic numbers, extra columns, and empty comment
    xyz = """
    3

    O1  0.0 0.0 0.0
    H1  1.0 0.0 0.0  0.5
    H2  0.0 1.0 0.0
    3
    second frame
    8   0.0 0.0 1.0
    1   1.0 0.0 1.0
    1   0.0 1.0 1.0

    """
    s = read_xyz(IOBuffer(xyz))
    @test length(s) == 6
    @test isnothing(s.title)
    @test name.(s) == ["O1", "H1", "H2", "8", "1", "1"]
    @test element.(s) == ["O", "H", "H", "O", "H", "H"]
    @test model.(s) == [1, 1, 1, 2, 2, 2]
    @test index.(s) == 1:6
    @test index_pdb.(s) == [1, 2, 3, 1, 2, 3]
    @test length(eachmodel(s)) == 2
    @test length(eachresidue(s)) == 2
    @test s[5].z == 1.0f0
    # errors
    @test_throws "expected the number of atoms" read_xyz(IOBuffer("abc\n"))
    @test_throws "expected 3 atoms" read_xyz(IOBuffer("3\n\nC 0 0 0\n"))
    @test_throws "could not read atom" read_xyz(IOBuffer("1\n\nC 0 0\n"))
    @test_throws "unknown atomic number" read_xyz(IOBuffer("1\n\n999 0 0 0\n"))
    @test_throws "Could not find any atom" read_xyz(IOBuffer("0\n\n"))
    @test_throws "Could not find any atom" read_xyz(PDBToolsBase.test_file("ethanol.xyz"), at -> false)
end
