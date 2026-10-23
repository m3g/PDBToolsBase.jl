"""
    read_mmcif(mmCIF_file::AbstractString, selection::AbstractString; field_assignment)
    read_mmcif(mmCIF_file::AbstractString, selection_function::Function = all, field_assignment)

    read_mmcif(mmCIF_data::IOBuffer, selection::AbstractString; field_assignment)
    read_mmcif(mmCIF_data::IOBuffer, selection_function::Function = all, field_assignment)

Reads a mmCIF file and stores the data in a vector of type `Atom`. 

All fields except the file name are optional.

If a selection is provided, only the atoms matching the selection will be read. 
For example, `resname ALA` will select all the atoms in the residue ALA.

If a selection function is provided, only the atoms for which `selection_function(atom)` is true will be returned.

The `field_assignment` keyword is `nothing` (default) or a `Dict{String,Symbol}` and can be used to specify which fields in the mmCIF file should be read into the `Atom` type.
For example `field_assignment = Dict("type_symbol" => :name)` will read the `_atom_site.type_symbol` field in the mmCIF 
file into the `name` field of the `Atom` type.

The default assignment is follows the standard mmCIF convention:

```julia
Dict{String,Symbol}(
    "id" => :index_pdb
    "Cartn_x" => :x
    "Cartn_y" => :y
    "Cartn_z" => :z
    "occupancy" => :occup
    "B_iso_or_equiv" => :beta
    "pdbx_formal_charge" => :charge
    "pdbx_PDB_model_num" => :model
    "label_atom_id" => :name
    "label_comp_id" => :resname
    "label_asym_id" => :chain
    "auth_seq_id" => :resnum
    "type_symbol" => :pdb_element
)
```

Source: https://mmcif.wwpdb.org/docs/tutorials/content/atomic-description.html

### Examples

```jldoctest
julia> using PDBTools

julia> ats = read_mmcif(PDBTools.TESTCIF)
   Vector{Atom{Nothing}} with 76 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     GLY     A        1        1   -4.564   25.503   24.113  1.00 24.33     1                 1
       2   CA     GLY     A        1        1   -4.990   26.813   24.706  1.00 24.29     1                 2
⋮
      75    O     HOH     Q       63       15   -3.585   34.725   20.903  1.00 19.82     1              2980
      76    O     HOH     Q       64       16   -4.799   40.689   37.419  1.00 20.13     1              2981

julia> ats = read_mmcif(PDBTools.TESTCIF, "index < 3")
   Vector{Atom{Nothing}} with 2 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     GLY     A        1        1   -4.564   25.503   24.113  1.00 24.33     1                 1
       2   CA     GLY     A        1        1   -4.990   26.813   24.706  1.00 24.29     1                 2

julia> ats = read_mmcif(PDBTools.TESTCIF, at -> name(at) == "CA")
   Vector{Atom{Nothing}} with 11 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       2   CA     GLY     A        1        1   -4.990   26.813   24.706  1.00 24.29     1                 2
       6   CA     GLN     A        2        2   -4.738   30.402   23.484  1.00 23.74     1                 6
⋮
      70   CA      CA     G     1003       10  -24.170   27.201   64.364  1.00 27.40     1              2967
      71   CA      CA     H     1004       11  -10.624   32.854   69.292  1.00 29.53     1              2968

```

"""
function read_mmcif end

function read_mmcif(cifdata::IOBuffer, selection_function::F=all; kargs...) where {F<:Function}
    atoms = _parse_mmCIF(cifdata, selection_function; kargs...)
    return atoms
end

function read_mmcif(filename::AbstractString, selection_function::F=all; kargs...) where {F<:Function}
    atoms = open(expanduser(filename), "r") do f
        _parse_mmCIF(f, selection_function; kargs...)
    end
    return atoms
end

function _maximum_read(atoms, stop_at, memory_available)
    if !isnothing(stop_at) && (length(atoms) >= stop_at)
        return true
    end
    if length(atoms) > 0 && mod(length(atoms), 50) == 0
        estimated_bytes = length(atoms) * (sizeof(eltype(atoms)) + sizeof(Ptr{Nothing}))
        if estimated_bytes > memory_available * Sys.total_memory()
            @warn """\n
                Memory limit reached. $(length(atoms)) atoms read so far will be returned.
                Size of the atoms array: $(round(estimated_bytes / 1024^2; digits=3)) MB

            """ _file = nothing _line = nothing
            return true
        end
    end
    return false
end

function _error_if_no_atoms(atoms)
    if length(atoms) == 0
        throw(ArgumentError("""\n 
            Could not find any atom in mmCIF file matching the selection. 

        """))
    end
end

# Internal function to replace the field assignment if a custom selection was
# provided by the user
function replace_custom_fields!(_atom_symbol_for_cif_field, field_assignment)
    if !isnothing(field_assignment)
        for (custom_key, custom_field) in field_assignment
            if !(custom_field in fieldnames(Atom))
                throw(ArgumentError("""\n
                    Field $custom_field not available in the PDBTools.Atom type.

                """))
            end
            if custom_field == :custom
                throw(ArgumentError("""\n
                    Setting a custom field to :custom is not currently supported.

                """))
            end
            # Type to assign to the field
            T = fieldtype(Atom, custom_field)
            # Remove the entry from the standard mmCIF fields, if available
            if haskey(_atom_symbol_for_cif_field, custom_key)
                pop!(_atom_symbol_for_cif_field, custom_key)
            end
            # Remove entry corresponding to custom_field, if available
            for p in _atom_symbol_for_cif_field
                key, val = p
                if last(val) == custom_field
                    pop!(_atom_symbol_for_cif_field, key)
                end
            end
            # Add entry associated to new assignment
            push!(_atom_symbol_for_cif_field, custom_key => (T, custom_field))
        end
    end
    return _atom_symbol_for_cif_field
end

#=
data_all
loop_
_atom_site.group_PDB             # ATOM 
_atom_site.id                    # 1 index_pdb
_atom_site.type_symbol           # N symbol?  
_atom_site.label_atom_id         # N name
_atom_site.label_alt_id          # . 
_atom_site.label_comp_id         # VAL resname
_atom_site.label_asym_id         # A chain
_atom_site.label_entity_id       # 1 ?
_atom_site.label_seq_id          # 1 resnum
_atom_site.pdbx_PDB_ins_code     # ? ? 
_atom_site.Cartn_x               # 6.204 x
_atom_site.Cartn_y               # 16.869 y
_atom_site.Cartn_z               # 4.854 z
_atom_site.occupancy             # 1.00 occup
_atom_site.B_iso_or_equiv        # 49.05 beta
_atom_site.pdbx_formal_charge    # ? charge
_atom_site.auth_seq_id           # 1
_atom_site.auth_comp_id          # VAL
_atom_site.auth_asym_id          # A
_atom_site.auth_atom_id          # N
_atom_site.pdbx_PDB_model_num    # 1 model
ATOM   1    N  N   . VAL A 1 1   ? 6.204   16.869  4.854   1.00 49.05 ? 1   VAL A N   1
ATOM   2    C  CA  . VAL A 1 1   ? 6.913   17.759  4.607   1.00 43.14 ? 1   VAL A CA  1
ATOM   3    C  C   . VAL A 1 1   ? 8.504   17.378  4.797   1.00 24.80 ? 1   VAL A C   1
ATOM   4    O  O   . VAL A 1 1   ? 8.805   17.011  5.943   1.00 37.68 ? 1   VAL A O   1
ATOM   5    C  CB  . VAL A 1 1   ? 6.369   19.044  5.810   1.00 72.12 ? 1   VAL A CB  1
ATOM   6    C  CG1 . VAL A 1 1   ? 7.009   20.127  5.418   1.00 61.79 ? 1   VAL A CG1 1
ATOM   7    C  CG2 . VAL A 1 1   ? 5.246   18.533  5.681   1.00 80.12 ? 1   VAL A CG2 1 
=#
function _supported_read_cif_fields(field_assignment)
    # We need this to be indexable (and ordered dict) to keep the order of the fields
    # when writing the mmCIF file
    _atom_symbol_for_cif_field = OrderedDict{String,Tuple{DataType,Symbol}}(
        "id" => (Int32, :index_pdb), # Standard mmCIF
        "type_symbol" => (StringType, :pdb_element), # Standard mmCIF
        "Cartn_x" => (Float32, :x),
        "Cartn_y" => (Float32, :y),
        "Cartn_z" => (Float32, :z),
        "occupancy" => (Float32, :occup),
        "B_iso_or_equiv" => (Float32, :beta),
        "pdbx_formal_charge" => (Float32, :charge),
        "auth_seq_id" => (Int32, :resnum), # Standard mmCIF
        "label_comp_id" => (StringType, :resname), # Standard mmCIF
        "label_asym_id" => (StringType, :chain), # Standard mmCIF
        "label_atom_id" => (StringType, :name), # Standard mmCIF
        "pdbx_PDB_model_num" => (Int32, :model),
    )
    return replace_custom_fields!(_atom_symbol_for_cif_field, field_assignment)
end

@testitem "_supported_read_cif_fields" setup = [PublicNames] begin
    using PDBToolsBase
    using PDBToolsBase: StringType
    using InlineStrings
    field_assignment = Dict(
        "type_symbol" => :name,
        "label_comp_id" => :resname,
        "label_asym_id" => :chain,
        "label_seq_id" => :resnum,
    )
    _cif_fields = PDBToolsBase._supported_read_cif_fields(field_assignment)
    @test length(_cif_fields) == 12
    @test _cif_fields["id"] == (Int32, :index_pdb)
    @test _cif_fields["type_symbol"] == (StringType, :name)
    @test _cif_fields["label_comp_id"] == (StringType, :resname)
    @test _cif_fields["label_asym_id"] == (StringType, :chain)
    @test _cif_fields["label_seq_id"] == (Int32, :resnum)
    @test _cif_fields["Cartn_x"] == (Float32, :x)
    @test _cif_fields["Cartn_y"] == (Float32, :y)
    @test _cif_fields["Cartn_z"] == (Float32, :z)
    @test _cif_fields["occupancy"] == (Float32, :occup)
    @test _cif_fields["B_iso_or_equiv"] == (Float32, :beta)
    @test _cif_fields["pdbx_formal_charge"] == (Float32, :charge)
    @test _cif_fields["pdbx_PDB_model_num"] == (Int32, :model)
    push!(field_assignment, "type_symbol" => :fail)
    @test_throws ArgumentError PDBToolsBase._supported_read_cif_fields(field_assignment)
    field_assignment = Dict("type_symbol" => :custom)
    @test_throws ArgumentError PDBToolsBase._supported_read_cif_fields(field_assignment)
    @test length(read_mmcif(PDBToolsBase.LONG_CHAIN_STRING_CIF)) == 8
    @test length(read_mmcif(PDBToolsBase.LONG_NAME_STRING_CIF)) == 1
    @test length(read_pdb(PDBToolsBase.LONG_CHAIN_STRING_CIF)) == 8
    @test length(read_pdb(PDBToolsBase.LONG_NAME_STRING_CIF)) == 1
end

function _new_cif_token(line)
    startswith(line, '_') && return true
    startswith(line, "loop_") && return true
    return false
end

function _parse_mmCIF(
    cifdata::Union{IOStream,IOBuffer}, selection_function::F;
    memory_available::Real=0.8,
    stop_at=nothing,
    field_assignment::Union{Nothing,Dict{String,Symbol}}=nothing,
) where {F<:Function}
    _atom_symbol_for_cif_field = _supported_read_cif_fields(field_assignment)
    atoms = Atom{Nothing}[]
    lastatom = Atom{Nothing}()
    # Columns of the _atom_site records to be read, and the Atom fields they are assigned to
    columns = Tuple{Int,Symbol}[]
    ncols = 0
    ifield = 0
    lattice = _empty_lattice
    _atom_site_field_inds = Dict{String,Int}()
    for line in eachline(cifdata)
        # Unit cell: the _cell fields are expected to precede the _atom_site fields
        if startswith(line, "_cell.")
            lattice = _read_lattice_line(lattice, line)
        end
        if startswith(line, "_atom_site.")
            field_end = findfirst(<=(' '), line)
            if isnothing(field_end)
                field_end = length(line) + 1
            end
            field = @view(line[12:field_end-1])
            ifield += 1
            _atom_site_field_inds[field] = ifield
        end
        # Header ended
        if ifield > 0 && startswith(line, r"ATOM|HETATM")
            for (key, keyval) in _atom_symbol_for_cif_field
                if haskey(_atom_site_field_inds, key)
                    push!(columns, (_atom_site_field_inds[key], last(keyval)))
                end
            end
            sort!(columns; by=first)
            ncols = length(keys(_atom_site_field_inds))
            if ncols > 0
                break
            end
        end
    end
    if length(columns) == 0
        throw(ArgumentError("Invalid mmCIF file: no ATOM or HETATM fields found."))
    end
    # Buffer for the fields of each record, reused for all atoms
    record_fields = SubString{String}[]
    sizehint!(record_fields, ncols)
    seekstart(cifdata)
    _atoms_token = false
    _atoms = false
    for line in eachline(cifdata)
        if startswith(line, "_atom_site.")
            _atoms_token = true
        end
        if _atoms_token 
            if startswith(line, r"ATOM|HETATM")
                _atoms = true
                atom = read_atom_mmcif(line, ncols, columns, record_fields, lastatom)
                selection_function(atom) && push!(atoms, atom)
                _maximum_read(atoms, stop_at, memory_available) && break
                lastatom = atom
            elseif _atoms # went through atom list and a new token is found
                _new_cif_token(line) && break
            end
        end
    end
    seekstart(cifdata)
    _error_if_no_atoms(atoms)
    return Structure(atoms, MetaData(_structure_unitcell(lattice)))
end

function read_atom_mmcif(
    record::String,
    ncols::Int,
    columns::Vector{Tuple{Int,Symbol}},
    record_fields::Vector{SubString{String}},
    lastatom::Atom,
)
    empty!(record_fields)
    for field in eachsplit(record)
        push!(record_fields, field)
    end
    if length(record_fields) != ncols
        throw(ArgumentError("mmCIF record has $(length(record_fields)) fields, expected $ncols: \"$record\""))
    end
    atom = Atom{Nothing}(; index=index(lastatom) + 1, residue=residue(lastatom))
    for (icol, field_name) in columns
        _set_atom_field!(atom, field_name, record_fields[icol])
    end
    if !same_residue(atom, lastatom)
        atom.residue = residue(lastatom) + 1
    end
    return atom
end

# Parses `value` into the field `field_name` of `atom`, with one branch per field of the
# Atom type, such that the type of each field is known at compile time (required for
# `juliac --trim`)
@generated function _set_atom_field!(atom::Atom, field_name::Symbol, value::AbstractString)
    ex = :(throw(ArgumentError(":$field_name is not a field of the Atom type that can be read from mmCIF files.")))
    for f in reverse(fieldnames(atom))
        f == :custom && continue
        T = fieldtype(atom, f)
        ex = :(field_name === $(QuoteNode(f)) ? setfield!(atom, $(QuoteNode(f)), _parse($T, value; alt=_alt($T))) : $ex)
    end
    return quote
        $ex
        return atom
    end
end

@testitem "read_mmcif" setup = [AllocTest, PublicNames] begin
    using PDBToolsBase
    using Chairmarks
    using .AllocTest: Allocs

    b = @b read_mmcif($(PDBToolsBase.TESTCIF))
    @test b.allocs < Allocs(500)
    ats = read_mmcif(PDBToolsBase.TESTCIF)
    @test count(iswater, ats) == 5
    @test count(isprotein, ats) == 69
    @test length(eachresidue(filter(isprotein, ats))) == 9
    tmpfile = tempname() * ".cif"
    write_mmcif(tmpfile, ats, isprotein)
    prot = read_mmcif(tmpfile)
    @test length(prot) == 69
    ats = read_mmcif(PDBToolsBase.TESTCIF, at -> resname(at) == "HOH")
    @test length(ats) == 5

    # performance tests for inner functions for reading cif fields into Atoms
    record = "ATOM   1    N  N   . VAL A 1 1   ? 6.204   16.869  4.854   1.00 49.05 ? 1   VAL A N   1"
    columns = [(2, :index_pdb), (4, :name), (6, :resname), (7, :chain), (9, :resnum), (11, :x), (12, :y), (13, :z), (14, :occup), (15, :beta), (16, :charge), (17, :resnum), (18, :resname), (19, :chain), (20, :name), (21, :model)]
    record_fields = SubString{String}[]
    lastatom = Atom()
    ncols = 21
    b = @b PDBToolsBase.read_atom_mmcif($record, $ncols, $columns, $record_fields, $lastatom)
    @test b.allocs == Allocs(1)
    atom = PDBToolsBase.read_atom_mmcif(record, ncols, columns, record_fields, lastatom)
    @test (name(atom), resname(atom), chain(atom), resnum(atom), model(atom)) == ("N", "VAL", "A", 1, 1)
    @test (atom.x, atom.y, atom.z, atom.occup, atom.beta, atom.charge) == (6.204f0, 16.869f0, 4.854f0, 1.0f0, 49.05f0, 0.0f0)
    b = @b PDBToolsBase._set_atom_field!($atom, :x, $(record_fields[11]))
    @test b.allocs == Allocs(0)
    @test_throws ArgumentError PDBToolsBase.read_atom_mmcif(record * " extra", ncols, columns, record_fields, lastatom)
    @test_throws ArgumentError PDBToolsBase._set_atom_field!(atom, :custom, record_fields[1])

    # Test early stoppers
    ats = read_mmcif(PDBToolsBase.TESTCIF; stop_at=10)
    @test length(ats) == 10
    ats = read_mmcif(PDBToolsBase.TESTCIF; memory_available=1e-10)
    @test length(ats) == 50

    # Read custom fields instead
    field_assignment = Dict(
        "type_symbol" => :name,
        "label_comp_id" => :resname,
        "label_asym_id" => :chain,
        "label_seq_id" => :resnum,
    )
    ats0 = read_mmcif(PDBToolsBase.TESTCIF)
    ats1 = read_mmcif(PDBToolsBase.TESTCIF; field_assignment)
    @test all(name.(ats1) .== pdb_element.(ats0))
    @test all(resname.(ats1) .== resname.(ats0))
    @test all(chain.(ats1) .== chain.(ats0))
    @test all(resnum.(filter(isprotein, ats1)) .== resnum.(filter(isprotein, ats0)))
    @test resnum.(filter(iswater, ats0)) == Int32[60, 61, 62, 63, 64]
    @test resnum.(filter(iswater, ats1)) == Int32[0, 0, 0, 0, 0]

    ats1 = read_pdb(PDBToolsBase.TESTCIF; field_assignment)
    @test all(name.(ats1) .== pdb_element.(ats0))
    @test all(resname.(ats1) .== resname.(ats0))
    @test all(chain.(ats1) .== chain.(ats0))
    @test all(resnum.(filter(isprotein, ats1)) .== resnum.(filter(isprotein, ats0)))
    @test resnum.(filter(iswater, ats0)) == Int32[60, 61, 62, 63, 64]
    @test resnum.(filter(iswater, ats1)) == Int32[0, 0, 0, 0, 0]

    # Throw error if no ATOM fields were found
    @test_throws "no ATOM or HETATM fields" read_mmcif(PDBToolsBase.BROKENCIF)
    @test_throws "no ATOM or HETATM fields" read_pdb(PDBToolsBase.BROKENCIF)
end