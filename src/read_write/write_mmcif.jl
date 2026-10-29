#=
loop_
_atom_site.group_PDB                  #  1 ATOM 
_atom_site.id                         #  2 index
_atom_site.type_symbol                #  3 element
_atom_site.label_atom_id              #  4 redundant: atom name (site.auth_atom_id)
_atom_site.label_alt_id               #  5 nothing
_atom_site.label_comp_id              #  6 redundant: residue name (site.auth_comp_id)
_atom_site.label_asym_id              #  7 redundant: chain id (site.auth_asym_id)
_atom_site.label_entity_id            #  8 entity_id: print 1 for protein, 2 for water, 3 for all others
_atom_site.label_seq_id               #  9 redundant: residue number (site.auth_seq_id)
_atom_site.pdbx_PDB_ins_code          # 10 nothing - print ?
_atom_site.Cartn_x                    # 11 x
_atom_site.Cartn_y                    # 12 y
_atom_site.Cartn_z                    # 13 z
_atom_site.occupancy                  # 14 occupancy
_atom_site.B_iso_or_equiv             # 15 bfactor
_atom_site.pdbx_formal_charge         # 16 charge: print ? if unknown
_atom_site.auth_seq_id                # 17 residue number   
_atom_site.auth_comp_id               # 18 residue name
_atom_site.auth_asym_id               # 19 chain id
_atom_site.auth_atom_id               # 20 atom name
_atom_site.pdbx_PDB_model_num         # 21 model number
ATOM   1    N  N   . VAL A 1 1   ? 6.204   16.869  4.854   1.00 49.05 ? 1   VAL A N   1
ATOM   2    C  CA  . VAL A 1 1   ? 6.913   17.759  4.607   1.00 43.14 ? 1   VAL A CA  1
ATOM   3    C  C   . VAL A 1 1   ? 8.504   17.378  4.797   1.00 24.80 ? 1   VAL A C   1
ATOM   4    O  O   . VAL A 1 1   ? 8.805   17.011  5.943   1.00 37.68 ? 1   VAL A O   1
ATOM   5    C  CB  . VAL A 1 1   ? 6.369   19.044  5.810   1.00 72.12 ? 1   VAL A CB  1
ATOM   6    C  CG1 . VAL A 1 1   ? 7.009   20.127  5.418   1.00 61.79 ? 1   VAL A CG1 1
ATOM   7    C  CG2 . VAL A 1 1   ? 5.246   18.533  5.681   1.00 80.12 ? 1   VAL A CG2 1 
=#
function _supported_write_cif_fields(field_assignment)
    # We need this to be indexable (and ordered dict) to keep the order of the fields
    # when writing the mmCIF file
    _atom_symbol_for_cif_field = OrderedDict{String,Tuple{DataType,Symbol}}(
        "id" => (Int32, :index_pdb), # Standard mmCIF
        "type_symbol" => (StringType, :pdb_element), # Standard mmCIF
        "label_atom_id" => (StringType, :name), # redundant - atom name
        "label_alt_id" => (String1, :_print_dot), # nothing (a dot: .)
        "label_comp_id" => (StringType, :resname), # redundant - residue name
        "label_asym_id" => (StringType, :chain), # redundant - chain id
        "label_entity_id" => (Int32, :_entity_id), # entity_id: print 1 for protein, 2 for water, 3 for all others
        "label_seq_id" => (Int32, :resnum), # redundant - residue number
        "pdbx_PDB_ins_code" => (String1, :_print_question_mark), # nothing - print ?
        "Cartn_x" => (Float32, :x),
        "Cartn_y" => (Float32, :y),
        "Cartn_z" => (Float32, :z),
        "occupancy" => (Float32, :occup),
        "B_iso_or_equiv" => (Float32, :beta),
        "pdbx_formal_charge" => (String1, :_print_question_mark), # use ? for unknown/zero charges
        "auth_seq_id" => (Int32, :resnum), # Standard mmCIF
        "auth_comp_id" => (StringType, :resname), # Standard mmCIF
        "auth_asym_id" => (StringType, :chain), # Standard mmCIF
        "auth_atom_id" => (StringType, :name), # Standard mmCIF
        "pdbx_PDB_model_num" => (Int32, :model),
    )
    return replace_custom_fields!(_atom_symbol_for_cif_field, field_assignment)
end

_write_mmcif_cell(::IO, ::Nothing) = nothing
function _write_mmcif_cell(io::IO, unitcell::AbstractMatrix)
    l = matrix_to_lattice(unitcell)
    println(io, @sprintf("_cell.length_a    %.3f", l.a))
    println(io, @sprintf("_cell.length_b    %.3f", l.b))
    println(io, @sprintf("_cell.length_c    %.3f", l.c))
    println(io, @sprintf("_cell.angle_alpha %.2f", l.α))
    println(io, @sprintf("_cell.angle_beta  %.2f", l.β))
    println(io, @sprintf("_cell.angle_gamma %.2f", l.γ))
    println(io, "#")
end

_dot_for_empty(val) = val
_dot_for_empty(val::AbstractString) = isnothing(findfirst(!isspace, val)) ? "." : val

# Formatting of the values written in the _atom_site loop. The kind (:int, :float or :string)
# is defined by the type of the field in the table of mmCIF fields.
const _fmt_string = Printf.Format("%$(_length(StringType))s ")
const _fmt_int = Printf.Format("%9i")
const _fmt_float = Printf.Format("%12.5f")
_fmt_kind(T::DataType) = T <: Integer ? :int : T <: AbstractFloat ? :float : :string
_fmt(::Symbol, val::AbstractString) = Printf.format(_fmt_string, val)
function _fmt(kind::Symbol, val::Real)
    kind === :int && return Printf.format(_fmt_int, val)
    kind === :float && return Printf.format(_fmt_float, val)
    return Printf.format(_fmt_string, val)
end

function _mmcif_field_string(atom::Atom, field_name::Symbol, kind::Symbol)
    field_name == :_print_atom && return _fmt(kind, "ATOM")
    field_name == :_print_dot && return _fmt(kind, ".")
    field_name == :_print_question_mark && return _fmt(kind, "?")
    field_name == :charge && return _fmt(kind, round(Int, charge(atom)))
    if field_name == :_entity_id
        isprotein(atom) && return _fmt(kind, 1)
        iswater(atom) && return _fmt(kind, 2)
        return _fmt(kind, 3)
    end
    return _format_atom_field(atom, field_name, kind)
end

# One branch per field of the Atom type, such that the type of each value is
# known at compile time (required for `juliac --trim`)
@generated function _format_atom_field(atom::Atom, field_name::Symbol, kind::Symbol)
    ex = :(throw(ArgumentError(":$field_name not found for printing mmcif field in atom")))
    for f in reverse(fieldnames(atom))
        f == :custom && continue
        ex = :(field_name === $(QuoteNode(f)) ? _fmt(kind, _dot_for_empty(getfield(atom, $(QuoteNode(f))))) : $ex)
    end
    return ex
end

"""
    write_mmcif(filename, atoms::AbstractVector{<:Atom}, [selection]; field_assignment=nothing, connectivity=true)

Write a mmCIF file with the atoms in `atoms` to `filename`. The optional `selection` argument is a string or function
that can be used to select a subset of the atoms in `atoms`. For example, `write_mmcif("test.cif", atoms, "name CA")`.

The optional `field_assignment` argument is a dictionary that can be used to assign custom fields to the mmCIF file.

If `connectivity=false`, the connectivity of the structure (`_geom_bond` fields) is not written.

If `atoms` is a `Structure`, its metadata is written to the file: the unit cell in the `_cell` 
fields, the title in the `_struct.title` field, and the connectivity in the `_geom_bond` 
fields, which list each bond once, by the `_atom_site.id` of the two atoms. The atoms are 
numbered sequentially in the file (`_atom_site.id` = 1, 2, ...), and the connectivity, which 
refers to the indices of the atoms (`index(atom)`), is converted to these numbers. Bonds to 
atoms that are not written (because of the selection) are not written. An error is thrown 
if two atoms to be written have the same index, because the bonds of these atoms would be
ambiguous: in this case, use `connectivity=false` to write the file without the connectivity.

"""
function write_mmcif(
    filename::AbstractString, 
    atoms::AbstractVector{<:Atom},
    selection_function::F=all;
    field_assignment::Union{Nothing,Dict{String,Symbol}}=nothing,
    connectivity::Bool=true,
) where {F<:Function}
    _cif_fields = _supported_write_cif_fields(field_assignment)
    bonds = connectivity ? _connectivity(atoms) : nothing
    # _atom_site.id of the atoms, by their indices (0 for atoms not written)
    atom_site_id = _written_serials(bonds, atoms, selection_function)
    open(expanduser(filename), "w") do file
        # Data block header
        println(file, "data_MODEL")
        println(file, "#")
        println(file, "_entry.id MODEL")
        println(file, "#")
        println(file, "_audit_conform.dict_name       mmcif_pdbx.dic")
        println(file, "_audit_conform.dict_version    5.397")
        println(file, "_audit_conform.dict_location   http://mmcif.pdb.org/dictionaries/ascii/mmcif_pdbx.dic")
        println(file, "#")
        _write_mmcif_title(file, _title(atoms))
        _write_mmcif_cell(file, _unitcell(atoms))
        # Collect chain→entity_id (unique per chain) and unique (chain, resnum) pairs
        chain_entity = OrderedDict{String,Int}()
        chain_type = OrderedDict{String,Symbol}()  # :protein, :water, :other
        residue_list = OrderedDict{Tuple{String,Int32},Nothing}()
        for atom in atoms
            !selection_function(atom) && continue
            ch = string(atom.chain)
            if !haskey(chain_entity, ch)
                chain_entity[ch] = length(chain_entity) + 1
                chain_type[ch] = isprotein(atom) ? :protein : iswater(atom) ? :water : :other
            end
            key = (ch, atom.resnum)
            haskey(residue_list, key) || (residue_list[key] = nothing)
        end
        # _entity: one entry per chain (each chain gets its own entity_id to avoid DSSP hang)
        println(file, "loop_")
        println(file, "_entity.id")
        println(file, "_entity.type")
        for (ch, eid) in chain_entity
            type_str = chain_type[ch] == :protein ? "polymer" : chain_type[ch] == :water ? "water" : "non-polymer"
            println(file, "$eid $type_str")
        end
        println(file, "#")
        # _entity_poly (only for polymer/protein chains)
        if any(t -> t == :protein, values(chain_type))
            println(file, "loop_")
            println(file, "_entity_poly.entity_id")
            println(file, "_entity_poly.type")
            for (ch, eid) in chain_entity
                chain_type[ch] == :protein && println(file, "$eid polypeptide(L)")
            end
            println(file, "#")
        end
        # _struct_asym
        println(file, "loop_")
        println(file, "_struct_asym.id")
        println(file, "_struct_asym.entity_id")
        for (ch, eid) in chain_entity
            println(file, "$ch $eid")
        end
        println(file, "#")
        # _pdbx_poly_seq_scheme (only for protein/polymer chains)
        poly_residues = [(ch, rnum) for (ch, rnum) in keys(residue_list) if get(chain_type, ch, :other) == :protein]
        if !isempty(poly_residues)
            println(file, "loop_")
            println(file, "_pdbx_poly_seq_scheme.asym_id")
            println(file, "_pdbx_poly_seq_scheme.seq_id")
            println(file, "_pdbx_poly_seq_scheme.pdb_strand_id")
            println(file, "_pdbx_poly_seq_scheme.pdb_seq_num")
            println(file, "_pdbx_poly_seq_scheme.pdb_ins_code")
            for (ch, rnum) in poly_residues
                println(file, "$ch $rnum $ch $rnum ?")
            end
            println(file, "#")
        end
        # _atom_site loop
        println(file, "loop_")
        println(file, "_atom_site.group_PDB")
        for (key, _) in _cif_fields
            println(file, "_atom_site.$key")
        end
        index = 0
        for atom in atoms
            !selection_function(atom) && continue
            index += 1
            buff = IOBuffer(; append=true)
            write(buff, "ATOM $(_fmt(:int, index))")
            for (_, field) in _cif_fields
                kind = _fmt_kind(first(field))
                field_name = last(field)
                field_name == :index_pdb && continue
                if field_name == :_entity_id
                    write(buff, _fmt(kind, chain_entity[string(atom.chain)]))
                else
                    write(buff, _mmcif_field_string(atom, field_name, kind))
                end
            end
            println(file, String(take!(buff)))
        end
        println(file, "#")
        _write_mmcif_geom_bond(file, bonds, atom_site_id)
    end
end

_write_mmcif_title(::IO, ::Nothing) = nothing
function _write_mmcif_title(io::IO, title::String)
    println(io, "_struct.entry_id MODEL")
    if !occursin('\'', title)
        println(io, "_struct.title '" * title * "'")
    elseif !occursin('"', title)
        println(io, "_struct.title \"" * title * "\"")
    else # text field
        println(io, "_struct.title")
        println(io, ";" * title)
        println(io, ";")
    end
    println(io, "#")
end

# Each bond is written once: the bond i-j is written for i < j, or for i > j if the 
# connectivity of j does not contain i
_write_mmcif_geom_bond(::IO, ::Nothing, atom_site_id) = nothing
function _write_mmcif_geom_bond(io::IO, connectivity::Vector{Vector{Int32}}, atom_site_id::Vector{Int32})
    header = false
    for (i, bonds) in pairs(connectivity)
        atom_site_id[i] == 0 && continue
        for j in bonds
            (1 <= j <= length(connectivity) && atom_site_id[j] != 0) || continue
            (i < j || !(i in connectivity[j])) || continue
            if !header
                println(io, "loop_")
                println(io, "_geom_bond.atom_site_id_1")
                println(io, "_geom_bond.atom_site_id_2")
                header = true
            end
            println(io, string(atom_site_id[i]) * " " * string(atom_site_id[j]))
        end
    end
    header && println(io, "#")
    return nothing
end

@testitem "write_mmcif" setup = [PublicNames] begin
    using PDBToolsBase
    ats = read_pdb(PDBToolsBase.SMALLPDB)
    tmpfile = tempname() * ".cif"
    write_mmcif(tmpfile, ats)
    @test isfile(tmpfile)
    ats_cif = read_mmcif(tmpfile)
    @test all(position(at1) ≈ position(at2) for (at1, at2) in zip(ats, ats_cif))
    field_assignment = Dict("test" => :name)
    write_mmcif(tmpfile, ats; field_assignment)
    ats0 = read_mmcif(tmpfile)
    @test all(name(at) == "X" for at in ats0)
    ats1 = read_mmcif(tmpfile; field_assignment)
    @test all(name(at1) == name(at2) for (at1, at2) in zip(ats, ats1))
    field_assignment = Dict("test" => :abc)
    @test_throws ArgumentError write_mmcif(tmpfile, ats; field_assignment)

    cA = filter(at -> name(at) == "CA", ats)
    write_mmcif(tmpfile, ats, at -> name(at) == "CA")
    @test isfile(tmpfile)
    ats_cif = read_mmcif(tmpfile)
    @test all(position(at1) ≈ position(at2) for (at1, at2) in zip(cA, ats_cif))

    at = Atom(name="HC", pdb_element="")
    write_mmcif(tmpfile, [at])
    ats = read_mmcif(tmpfile)
    @test element(ats[1]) == "H"

    # unit cell of Structures is written in the _cell fields
    uc = lattice_to_matrix(54.53, 63.78, 82.34, 72.01, 89.44, 89.43)
    ats = read_pdb(PDBToolsBase.SMALLPDB)
    write_mmcif(tmpfile, Structure(ats; unitcell=uc))
    @test read_unitcell(tmpfile) ≈ uc rtol = 1e-4
    @test all(position(at1) ≈ position(at2) for (at1, at2) in zip(ats, read_mmcif(tmpfile)))
    write_mmcif(tmpfile, Structure(ats))
    @test_throws ArgumentError read_unitcell(tmpfile)

end

@testitem "write_mmcif - connectivity and title" setup = [PublicNames] begin
    using PDBToolsBase
    file = joinpath(pkgdir(PDBToolsBase), "test", "ethanol.pdb")
    s = read_pdb(file)
    tmpfile = tempname() * ".cif"
    write_mmcif(tmpfile, s)
    lines = readlines(tmpfile)
    @test "_struct.title 'ETHANOL MOLECULE WITH EXPLICIT HYDROGENS'" in lines
    ibond = findfirst(==("_geom_bond.atom_site_id_2"), lines)
    @test lines[ibond+1:ibond+9] == ["1 2", "1 4", "1 5", "1 6", "2 3", "2 7", "2 8", "3 9", "#"]
    s2 = read_mmcif(tmpfile)
    @test s2.title == s.title
    @test s2.connectivity == s.connectivity
    @test read_pdb(tmpfile).connectivity == s.connectivity
    # selection: the atoms are renumbered in the file
    write_mmcif(tmpfile, s, at -> name(at) in ("C2", "O1", "H6"))
    lines = readlines(tmpfile)
    ibond = findfirst(==("_geom_bond.atom_site_id_2"), lines)
    @test lines[ibond+1:ibond+3] == ["1 2", "2 3", "#"]
    s2 = read_mmcif(tmpfile)
    @test s2.connectivity == [[2], [1, 3], [2]]
    # bonds listed in one direction only
    s.connectivity = [[2], Int32[], [2]]
    write_mmcif(tmpfile, s)
    lines = readlines(tmpfile)
    ibond = findfirst(==("_geom_bond.atom_site_id_2"), lines)
    @test lines[ibond+1:ibond+3] == ["1 2", "3 2", "#"]
    @test read_mmcif(tmpfile).connectivity[1:3] == [[2], [1, 3], [2]]
    # titles with quotes
    for title in ("it's a title", "it's a \"title\"", "a \"title\"")
        write_mmcif(tmpfile, Structure(s.atoms; title))
        @test read_mmcif(tmpfile).title == title
    end
    # no metadata
    write_mmcif(tmpfile, s.atoms)
    s2 = read_mmcif(tmpfile)
    @test isnothing(s2.title) && isnothing(s2.connectivity)
end
