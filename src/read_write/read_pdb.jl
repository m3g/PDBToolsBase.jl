"""
    read_pdb(pdbfile::AbstractString, selection::AbstractString)
    read_pdb(pdbfile::AbstractString, selection_function::Function = all)

    read_pdb(pdbdata::IOBuffer, selection::AbstractString)
    read_pdb(pdbdata::IOBuffer, selection_function::Function = all)

Reads a PDB file (in legacy PDB or mmCIF/PDBx formats) and stores the data in a vector of type `Atom`. 

If a selection is provided, only the atoms matching the selection will be read. 
For example, `resname ALA` will select all the atoms in the residue ALA.

If a selection function keyword is provided, only the atoms for which `selection_function(atom)` is true will be read.

The function automatically detects mmCIF format and switches to the mmCIF reader if the `loop_` keyword is found in the file, 
or if an `ATOM`/`HETATM` line cannot be parsed in PDB format.

When reading an mmCIF file, the following keyword arguments are forwarded to the mmCIF reader:
- `field_assignment::Union{Nothing,Dict{String,Symbol}} = nothing`: custom mapping from mmCIF `_atom_site` field names to `Atom` fields. See [`read_mmcif`](@ref) for details.
- `memory_available::Real = 0.8`: fraction of available memory the reader is allowed to use.
- `stop_at = nothing`: stop reading after this many atoms (useful for large files).

### Examples

```jldoctest
julia> using PDBTools

julia> protein = read_pdb(PDBTools.test_file("structure.pdb"))
   Structure{Atom{Nothing}} with 62026 atoms
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
⋮
   62025   H1    TIP3     C     9339    19638   13.218   -3.647  -34.453  1.00  0.00     1    WAT2     62025
   62026   H2    TIP3     C     9339    19638   12.618   -4.977  -34.303  1.00  0.00     1    WAT2     62026

julia> ALA = read_pdb(PDBTools.test_file("structure.pdb"),"resname ALA")
   Structure{Atom{Nothing}} with 72 atoms
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
⋮
    1339    C     ALA     A       95       95   14.815   -3.057   -5.633  1.00  0.00     1    PROT      1339
    1340    O     ALA     A       95       95   14.862   -2.204   -6.518  1.00  0.00     1    PROT      1340

julia> ALA = read_pdb(PDBTools.test_file("structure.pdb"), atom -> atom.resname == "ALA")
   Structure{Atom{Nothing}} with 72 atoms
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  0.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
⋮
    1339    C     ALA     A       95       95   14.815   -3.057   -5.633  1.00  0.00     1    PROT      1339
    1340    O     ALA     A       95       95   14.862   -2.204   -6.518  1.00  0.00     1    PROT      1340
```

"""
function read_pdb end

function _switch_to_mmcif(pdbdata::Union{IOStream, IOBuffer})
    switch = false
    for line in eachline(pdbdata)
        if startswith(line, "loop_")
            switch = true
            break
        else
            try
                if startswith(line, r"ATOM|HETATM")
                    read_atom_pdb(line)
                    break
                end
            catch
                @warn "Failed reading ATOM line in PDB format, trying mmCIF reader." _line = nothing _file = nothing
                switch = true
                break
            end
        end
    end
    seekstart(pdbdata)
    return switch
end

function read_pdb(pdbdata::IOBuffer, selection_function::F=all; kargs...) where {F<:Function}
    atoms = if _switch_to_mmcif(pdbdata) 
        _parse_mmCIF(pdbdata, selection_function; kargs...)
    else
        _parse_pdb(pdbdata, selection_function)
    end
    return atoms
end

function read_pdb(filename::AbstractString, selection_function::F=all; kargs...) where {F<:Function}
    atoms = open(expanduser(filename), "r") do f
        if _switch_to_mmcif(f)
            _parse_mmCIF(f, selection_function; kargs...)
        else
            _parse_pdb(f, selection_function)
        end
    end
    return atoms
end

function _parse_pdb(pdbdata::Union{IOStream,IOBuffer}, selection_function::F) where {F<:Function}
    imodel = 1
    atoms = Atom{Nothing}[]
    lastatom = Atom{Nothing}()
    lattice = _empty_lattice
    for line in eachline(pdbdata)
        if occursin("END", line)
            imodel = imodel + 1
        end
        # Unit cell: from the first CRYST1 field
        if startswith(line, "CRYST1") && !_is_complete(lattice)
            lattice = _read_lattice_line(lattice, line)
        end
        if startswith(line, r"ATOM|HETATM")
            atom = read_atom_pdb(line, lastatom, imodel)
            atom.model = imodel
            selection_function(atom) && push!(atoms, atom)
            lastatom = atom
        end
    end
    seekstart(pdbdata)
    if length(atoms) == 0
        throw(ArgumentError("""\n 
            Could not find any atom in PDB file matching the selection. 

        """))
    end
    return Structure(atoms, MetaData(_structure_unitcell(lattice)))
end

# read atom from PDB file
function read_atom_pdb(record::AbstractString, lastatom::Atom=Atom{Nothing}(), imodel::Int=1)
    atom = Atom{Nothing}(; index=index(lastatom) + 1, residue=residue(lastatom), model=imodel, flag=0)
    record[1:6] == "HETATM" && (atom.flag = 1)
    inds_and_names = (
        (1, Val(:name)),
        (2, Val(:resname)),
        (3, Val(:chain)),
        (4, Val(:index_pdb)),
        (5, Val(:resnum)),
        (6, Val(:x)),
        (7, Val(:y)),
        (8, Val(:z)),
        (9, Val(:occup)),
        (10, Val(:beta)),
        (11, Val(:segname)),
        (12, Val(:pdb_element)),
        (13, Val(:charge)),
    )
    @views field_values = (
        record[13:16], # :name
        record[17:21], # :resname
        record[22:22], # :chain
        record[7:11],  # :index_pdb
        record[23:26], # :resnum
        record[31:38], # :x
        record[39:46], # :y
        record[47:54], # :z
        length(record) >= 60 ? record[56:60] : record[56:end], # :occup
        length(record) >= 66 ? record[61:66] : record[61:end], # :beta
        length(record) >= 76 ? record[73:76] : record[73:end], # :segname
        length(record) >= 78 ? record[77:78] : record[77:end], # :pdb_element
        length(record) >= 80 ? _parse_charge(record[79:80]) : _parse_charge(record[79:end]), # :charge 
    )
    _fast_setfield!(atom, field_values, inds_and_names)
    if !same_residue(atom, lastatom)
        atom.residue = residue(lastatom) + 1
    end
    return atom
end

@testitem "read_pdb" setup = [PublicNames] begin
    pdb_file = PDBToolsBase.test_file("structure.pdb")
    pdbdata = read(pdb_file, String)
    atoms = read_pdb(IOBuffer(pdbdata), at -> isprotein(at) && name(at) == "CA")
    @test length(atoms) == 104
    atoms = read_pdb(pdb_file, at -> isprotein(at) && name(at) == "CA")
    @test length(atoms) == 104

    cif_file = PDBToolsBase.TESTCIF
    pdbdata = read(cif_file, String)
    atoms = read_pdb(IOBuffer(pdbdata), at -> isprotein(at) && name(at) == "CA")
    atoms_ref = read_mmcif(cif_file, at -> isprotein(at) && name(at) == "CA")
    @test atoms == atoms_ref
end

@testitem "read_atom_pdb" setup = [PublicNames] begin
    line = "HETATM    1  O1  BGL    1       0.665   1.214  -0.259  1.00  0.00"
    a = PDBToolsBase.read_atom_pdb(line)
    @test a.index == 1
    @test a.name == "O1"
    @test a.resname == "BGL"
    @test a.chain == "X"
    @test a.resnum == 1
    @test a.residue == 1
    @test (a.x, a.y, a.z) == (0.665f0, 1.214f0, -0.259f0)
    @test (a.occup, a.beta) == (1.0, 0.0)
    @test a.model == 1
    @test a.segname == "X"
    @test a.index_pdb == 1
    @test a.pdb_element == "X"
    @test a.charge == 0.0f0

    line = "HETATM    1  O1  AGL     1       0.803   1.186  -0.211  1.00  0.00            1+  "
    a = PDBToolsBase.read_atom_pdb(line)
    @test a.index == 1
    @test a.name == "O1"
    @test a.resname == "AGL"
    @test a.chain == "X"
    @test a.resnum == 1
    @test a.residue == 1
    @test (a.x, a.y, a.z) == (0.803f0, 1.186f0, -0.211f0)
    @test (a.occup, a.beta) == (1.0, 0.0)
    @test a.model == 1
    @test a.segname == "X"
    @test a.index_pdb == 1
    @test a.pdb_element == "X"
    @test a.charge == 1.0f0

    line = "HETATM    1  O1  AGL     1       0.803   1.186  -0.211  1.00  0.00            +1  "
    a = PDBToolsBase.read_atom_pdb(line)
    @test a.charge == 1.0f0

    line = "HETATM    1  O1  AGL     1       0.803   1.186  -0.211  1.00  0.00            -1  "
    a = PDBToolsBase.read_atom_pdb(line)
    @test a.charge == -1.0f0

    line = "HETATM    1  O1  AGL     1       0.803   1.186  -0.211  1.00  0.00            1-  "
    a = PDBToolsBase.read_atom_pdb(line)
    @test a.charge == -1.0f0

    line = "ATOM  *****  H2  GLYCD4301      11.014  36.823  40.115  1.00  0.00      GLYC H"
    a = PDBToolsBase.read_atom_pdb(line)
    @test a.index == 1
    @test a.name == "H2"
    @test a.resname == "GLYC"
    @test a.chain == "D"
    @test a.resnum == 4301
    @test a.residue == 1
    @test (a.x, a.y, a.z) == (11.014f0, 36.823f0, 40.115f0)
    @test (a.occup, a.beta) == (1.0, 0.0)
    @test a.model == 1
    @test a.segname == "GLYC"
    @test a.index_pdb == 0
    @test a.pdb_element == "H"
    @test a.charge == 0

    line = "ATOM  3ce2c  LP2 WLS Cffff     376.512 638.670  16.990  0.00  0.00         C"
    a = PDBToolsBase.read_atom_pdb(line)
    @test a.index == 1
    @test a.index_pdb == 249388
    @test a.resnum == 65535
end

@testitem "abstract string" setup = [PublicNames] begin
    file = PDBToolsBase.test_file("structure.pdb") * "abc"
    s = read_pdb(@view(file[1:end-3]))
    @test length(s) == 62026
    tmpfile = tempname()
    write_pdb(tmpfile, s)
    @test isfile(tmpfile)
end

@testitem "switch to mmCIF" setup = [PublicNames] begin
    p1 = read_mmcif(PDBToolsBase.test_file("6CO8.cif"))
    p2 = read_pdb(PDBToolsBase.test_file("6CO8.cif"))
    @test p1 == p2

    p1 = read_mmcif(PDBToolsBase.TESTCIF)
    p2 = read_pdb(PDBToolsBase.TESTCIF)
    @test p1 == p2
end

@testitem "read Structure and unit cell" setup = [PublicNames] begin
    using PDBToolsBase
    # PDB file with CRYST1 field
    s = read_pdb(PDBToolsBase.test_file("pbc_protein.pdb"))
    @test s isa Structure{Atom{Nothing}}
    @test s.unitcell isa UnitCell
    @test s.unitcell ≈ read_unitcell(PDBToolsBase.test_file("pbc_protein.pdb"))
    @test read_pdb(IOBuffer(read(PDBToolsBase.test_file("pbc_protein.pdb"), String))).unitcell == s.unitcell
    @test read_pdb(PDBToolsBase.test_file("pbc_protein.pdb"), isprotein).unitcell == s.unitcell
    # mmCIF files with _cell fields
    s = read_mmcif(PDBToolsBase.TESTCIF)
    @test s isa Structure{Atom{Nothing}}
    @test s.unitcell ≈ read_unitcell(PDBToolsBase.TESTCIF)
    @test read_pdb(PDBToolsBase.TESTCIF).unitcell == s.unitcell
    @test read_mmcif(IOBuffer(read(PDBToolsBase.TESTCIF, String))).unitcell == s.unitcell
    lbd = PDBToolsBase.test_file("1LBD.cif")
    @test matrix_to_lattice(read_mmcif(lbd).unitcell).γ ≈ 120
    # files without unit cell
    @test isnothing(read_pdb(PDBToolsBase.SMALLPDB).unitcell)
    @test isnothing(read_mmcif(PDBToolsBase.test_file("6CO8.cif")).unitcell)
    # the 1Å cube indicates a non-crystallographic structure (here, cryo-EM)
    em = PDBToolsBase.test_file("7YWH.cif")
    @test read_unitcell(em) ≈ lattice_to_matrix(1.0, 1.0, 1.0, 90.0, 90.0, 90.0)
    @test isnothing(read_mmcif(em).unitcell)
    # round trip
    tmp = tempname()
    s = read_pdb(PDBToolsBase.test_file("pbc_protein.pdb"))
    write_pdb(tmp * ".pdb", s)
    @test read_pdb(tmp * ".pdb").unitcell ≈ s.unitcell
    write_mmcif(tmp * ".cif", s)
    @test read_mmcif(tmp * ".cif").unitcell ≈ s.unitcell
end
