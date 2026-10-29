"""
    write_pdb(filename::AbstractString, atoms::AbstractVector{<:Atom}, [selection]; header=:auto, footer=:auto, append=false, connectivity=true)

Write a PDB file with the atoms in `atoms` to `filename`. The `selection` argument is a string or function
that can be used to select a subset of the atoms in `atoms`. For example, `write_pdb("test.pdb", atoms, "name CA")`.

# Arguments

- `filename::AbstractString`: The name of the file to write.
- `atoms::AbstractVector{<:Atom}`: The atoms to write to the file.

# Optional positional argument 

- `selection::AbstractString`: A selection string to select a subset of the atoms in `atoms`.

# Keyword arguments

- `header::Union{String, Nothing}=:auto`: The header to add to the PDB file. If `:auto`, a header will be added with the number of atoms in `atoms`.
- `footer::Union{String, Nothing}=:auto`: The footer to add to the PDB file. If `:auto`, a footer will be added with the "END" keyword.
- `append::Bool=false`: If `true`, the atoms will be appended to the file instead of overwriting it.
- `connectivity::Bool=true`: If `false`, the connectivity of the structure (`CONECT` records) is not written.

The atoms are numbered sequentially in the file (serial numbers 1, 2, ..., in hexadecimal
format above 99999), in the order in which they are written.

If `atoms` is a `Structure`, its metadata is written to the file: the unit cell in the `CRYST1`
field, the title in `TITLE` records, and the connectivity in `CONECT` records. The 
connectivity, which refers to the indices of the atoms (`index(atom)`), is converted to the 
serial numbers of the atoms in the file. Bonds to atoms that are not written (because of the 
selection) are not written. An error is thrown if two atoms to be written have the same index,
because the bonds of these atoms would be ambiguous: in this case, use `connectivity=false` 
to write the file without the connectivity.

"""
function write_pdb(
    filename::AbstractString, 
    atoms::AbstractVector{<:Atom},
    selection_function::F=all;
    append=false,
    header=:auto,
    footer=:auto,
    connectivity::Bool=true,
) where {F<:Function}
    bonds = connectivity ? _connectivity(atoms) : nothing
    serials = _written_serials(bonds, atoms, selection_function)
    open(expanduser(filename), append ? "a" : "w") do io
        _write_pdb_header(io, header, length(atoms))
        _write_pdb_title(io, _title(atoms))
        _write_pdb_cryst1(io, _unitcell(atoms))
        serial = 0
        for atom in atoms
            if selection_function(atom)
                serial += 1
                println(io, write_pdb_atom(atom, serial))
            end
        end
        _write_pdb_conect(io, bonds, serials, atoms, selection_function)
        _write_pdb_footer(io, footer)
    end
end

# Header and footer writers dispatch on the type of the keyword values, which
# keeps `write_pdb` type-stable (required for `juliac --trim`)
const _pdb_header_date_format = Dates.DateFormat("dd-u-yy")
_write_pdb_header(::IO, ::Nothing, natoms::Int) = nothing
_write_pdb_header(io::IO, header, natoms::Int) = println(io, header)
function _write_pdb_header(io::IO, header::Symbol, natoms::Int)
    if header === :auto
        curr_date = Dates.format(Dates.today(), _pdb_header_date_format)
        println(io, @sprintf "%-10s%-40s%9s" "HEADER" "PDBTools.jl - $natoms atoms" curr_date)
    else
        println(io, header)
    end
end
_write_pdb_cryst1(::IO, ::Nothing) = nothing
function _write_pdb_cryst1(io::IO, unitcell::AbstractMatrix)
    l = matrix_to_lattice(unitcell)
    println(io, @sprintf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4i", l.a, l.b, l.c, l.α, l.β, l.γ, "P 1", 1))
end

# TITLE records: the title is written in columns 11-80, with continuation lines numbered 
# in columns 9-10. Lines are broken at spaces, when possible.
_write_pdb_title(::IO, ::Nothing) = nothing
function _write_pdb_title(io::IO, title::String)
    words = split(title)
    isempty(words) && return nothing
    iline = 1
    line = ""
    for word in words
        width = iline == 1 ? 70 : 69 # continuation lines start with a space
        if !isempty(line) && length(line) + 1 + length(word) > width
            _print_pdb_title_line(io, iline, line)
            iline += 1
            line = ""
        end
        # words longer than the line are split
        while length(word) > (iline == 1 ? 70 : 69)
            n = iline == 1 ? 70 : 69
            _print_pdb_title_line(io, iline, first(word, n))
            iline += 1
            word = @view(word[nextind(word, 0, n + 1):end])
        end
        line = isempty(line) ? String(word) : line * " " * word
    end
    _print_pdb_title_line(io, iline, line)
    return nothing
end
function _print_pdb_title_line(io::IO, iline::Int, text::AbstractString)
    if iline == 1
        println(io, "TITLE     " * text)
    else
        println(io, @sprintf("TITLE   %2i %s", iline, text))
    end
end

# Serial numbers of the atoms in the written file (1, 2, ..., in the order in which the atoms
# are written), by their indices: `serials[index(atom)]` is the serial number of the atom, or 
# 0 if the atom is not written. Only the indices covered by the connectivity are mapped.
_written_serials(::Nothing, atoms, selection_function) = Int32[]
function _written_serials(
    connectivity::Vector{Vector{Int32}}, atoms::AbstractVector{<:Atom}, selection_function::F,
) where {F<:Function}
    serials = zeros(Int32, length(connectivity))
    serial = 0
    for atom in atoms
        selection_function(atom) || continue
        serial += 1
        i = index(atom)
        1 <= i <= length(serials) || continue
        if serials[i] != 0
            throw(ArgumentError("""\n
                Cannot write the connectivity: more than one atom has index $i. 
                The connectivity refers to the indices of the atoms (`index(atom)`), which must be unique.
                Renumber the atoms (and their connectivity), or use `connectivity=false` to write the
                file without the connectivity.

            """))
        end
        serials[i] = serial
    end
    return serials
end

# CONECT records: up to 4 bonded atoms per record, with the serial numbers written as in
# the ATOM records (hexadecimal above 99999)
_write_pdb_conect(::IO, ::Nothing, serials, atoms, selection_function) = nothing
function _write_pdb_conect(
    io::IO, connectivity::Vector{Vector{Int32}}, serials::Vector{Int32}, 
    atoms::AbstractVector{<:Atom}, selection_function::F,
) where {F<:Function}
    bonded = Int32[]
    for atom in atoms
        selection_function(atom) || continue
        i = index(atom)
        1 <= i <= length(connectivity) || continue
        empty!(bonded)
        for j in connectivity[i]
            1 <= j <= length(serials) && serials[j] != 0 && push!(bonded, serials[j])
        end
        for ifirst in 1:4:length(bonded)
            print(io, "CONECT", _pdb_serial(serials[i]))
            for j in @view(bonded[ifirst:min(ifirst + 3, length(bonded))])
                print(io, _pdb_serial(j))
            end
            println(io)
        end
    end
    return nothing
end
_pdb_serial(i::Integer) = i <= 99999 ? @sprintf("%5i", i) : @sprintf("%5x", mod(i, 1048576))

_write_pdb_footer(::IO, ::Nothing) = nothing
_write_pdb_footer(io::IO, footer) = println(io, footer)
_write_pdb_footer(io::IO, footer::Symbol) = println(io, footer === :auto ? "END" : String(footer))

@testitem "write_pdb" setup = [PublicNames] begin
    using PDBToolsBase
    using DelimitedFiles
    pdb = read_pdb(PDBToolsBase.SMALLPDB)
    tmpfile = tempname() * ".pdb"
    write_pdb(tmpfile, pdb)
    @test isfile(tmpfile)
    f1 = readdlm(PDBToolsBase.SMALLPDB, '\n', header=true)
    f2 = readdlm(tmpfile, '\n', header=true)
    @test f1[1] == f2[1]
    # test selection
    write_pdb(tmpfile, pdb, at -> name(at) == "CA")
    f3 = read_pdb(tmpfile)
    @test length(f3) == 3
    # test header and footer
    write_pdb(tmpfile, pdb, at -> name(at) == "CA"; header="HEADER test", footer="END test")
    s = split(String(read(tmpfile)))
    @test (s[begin], s[begin+1]) == ("HEADER", "test")
    @test (s[end-1], s[end]) == ("END", "test")
    # test append
    append_pdb = tempname() * ".pdb"
    write_pdb(append_pdb, pdb; append=true)
    pdb1 = read_pdb(append_pdb)
    @test length(eachmodel(pdb1)) == 1
    @test length(pdb1) == 35
    write_pdb(append_pdb, pdb; append=true)
    pdb2 = read_pdb(append_pdb)
    @test length(eachmodel(pdb2)) == 2
    @test length(pdb2) == 70
    pdb_with_hetatm = read_pdb(PDBToolsBase.test_file("hetatm.pdb"))
    @test length(pdb_with_hetatm) == 19
    @test length(eachresidue(pdb_with_hetatm)) == 3
    @test getfield.(pdb_with_hetatm, :flag) == [i < 12 ? 0 : 1 for i in 1:length(pdb_with_hetatm)]
    tmpfile = tempname() * ".pdb"
    write_pdb(tmpfile, pdb_with_hetatm)
    pdb_with_hetatm = read_pdb(tmpfile)
    @test length(pdb_with_hetatm) == 19
    @test length(eachresidue(pdb_with_hetatm)) == 3
    @test getfield.(pdb_with_hetatm, :flag) == [i < 12 ? 0 : 1 for i in 1:length(pdb_with_hetatm)]

    pdb = read_pdb(PDBToolsBase.SMALLPDB)
    tmpfile = tempname() * ".pdb"
    write_pdb(tmpfile, pdb, at -> name(at) == "CA")
    @test isfile(tmpfile)
    ats = read_pdb(tmpfile)
    @test length(ats) == 3

    # unit cell of Structures is written in the CRYST1 field
    uc = lattice_to_matrix(54.53, 63.78, 82.34, 72.01, 89.44, 89.43)
    write_pdb(tmpfile, Structure(pdb; unitcell=uc))
    @test count(startswith("CRYST1"), readlines(tmpfile)) == 1
    @test read_unitcell(tmpfile) ≈ uc rtol = 1e-4
    @test length(read_pdb(tmpfile)) == length(pdb)
    write_pdb(tmpfile, Structure(pdb); header=nothing)
    @test !any(startswith("CRYST1"), readlines(tmpfile))
end

@testitem "write_pdb - connectivity and title" setup = [PublicNames] begin
    using PDBToolsBase
    file = PDBToolsBase.test_file("ethanol.pdb")
    s = read_pdb(file)
    tmpfile = tempname() * ".pdb"
    write_pdb(tmpfile, s)
    lines = readlines(tmpfile)
    @test lines[2] == "TITLE     ETHANOL MOLECULE WITH EXPLICIT HYDROGENS"
    @test count(startswith("CONECT"), lines) == 9
    @test filter(startswith("CONECT"), lines) == filter(startswith("CONECT"), readlines(file))
    @test last(lines) == "END"
    s2 = read_pdb(tmpfile)
    @test s2.title == s.title
    @test s2.connectivity == s.connectivity
    # selection: bonds to atoms not written are not written
    write_pdb(tmpfile, s, at -> element(at) != "H")
    @test filter(startswith("CONECT"), readlines(tmpfile)) == ["CONECT    1    2", "CONECT    2    1    3", "CONECT    3    2"]
    s2 = read_pdb(tmpfile)
    @test s2.connectivity == [[2], [1, 3], [2]]
    # the atoms are renumbered in the file, and the connectivity is converted
    write_pdb(tmpfile, s, at -> name(at) in ("C2", "O1", "H6"))
    lines = readlines(tmpfile)
    @test [l[7:11] for l in filter(startswith("HETATM"), lines)] == ["    1", "    2", "    3"]
    @test filter(startswith("CONECT"), lines) == ["CONECT    1    2", "CONECT    2    1    3", "CONECT    3    2"]
    # bonds are preserved when the order of the atoms changes
    bonded_names(s) = [sort([name(s[findfirst(at -> index(at) == j, s)]) for j in s.connectivity[index(at)]]) for at in s]
    r = reverse(s)
    write_pdb(tmpfile, r)
    r2 = read_pdb(tmpfile)
    @test name.(r2) == name.(r)
    @test bonded_names(r2) == bonded_names(r)
    @test r2.connectivity[1] == [7] # H6 - O1
    # repeated indices
    @test_throws "more than one atom has index" write_pdb(tmpfile, vcat(s, s))
    @test_throws "more than one atom has index" write_mmcif(tmpfile, vcat(s, s))
    @test_throws "connectivity=false" write_pdb(tmpfile, vcat(s, s))
    # skipping the connectivity
    write_pdb(tmpfile, vcat(s, s); connectivity=false)
    s2 = read_pdb(tmpfile)
    @test length(s2) == 18
    @test isnothing(s2.connectivity)
    @test s2.title == s.title
    write_mmcif(tmpfile, vcat(s, s); connectivity=false)
    s2 = read_mmcif(tmpfile)
    @test length(s2) == 18
    @test isnothing(s2.connectivity)
    @test s2.title == s.title
    write_pdb(tmpfile, vcat(s, s).atoms) # no connectivity: no error
    @test length(read_pdb(tmpfile)) == 18
    # more than 4 bonds, and indices greater than 99999
    at = Atom(name="C")
    s = Structure([at]; connectivity=[[2, 3, 4, 5, 6], [1], [1], [1], [1], [1]])
    io = IOBuffer()
    ats = [Atom(index=i) for i in 1:6]
    serials = PDBToolsBase._written_serials(s.connectivity, ats, PDBToolsBase.all)
    PDBToolsBase._write_pdb_conect(io, s.connectivity, serials, ats, PDBToolsBase.all)
    @test String(take!(io)) == "CONECT    1    2    3    4    5\nCONECT    1    6\n" * join(["CONECT    $i    1\n" for i in 2:6])
    @test PDBToolsBase._pdb_serial(100000) == "186a0"
    # long titles
    title = join(fill("ABCDEFGHI", 20), " ")
    write_pdb(tmpfile, Structure(read_pdb(file); title))
    lines = filter(startswith("TITLE"), readlines(tmpfile))
    @test length(lines) == 3
    @test all(length(l) <= 80 for l in lines)
    @test startswith(lines[2], "TITLE    2 ABCDEFGHI")
    @test read_pdb(tmpfile).title == title
    title = repeat("A", 150)
    write_pdb(tmpfile, Structure(read_pdb(file); title))
    lines = filter(startswith("TITLE"), readlines(tmpfile))
    @test length(lines) == 3 && all(length(l) <= 80 for l in lines)
    # vectors of atoms (not structures), and structures without title and connectivity
    write_pdb(tmpfile, read_pdb(file).atoms)
    @test !any(l -> startswith(l, "CONECT") || startswith(l, "TITLE"), readlines(tmpfile))
    write_pdb(tmpfile, Structure(read_pdb(file).atoms))
    @test !any(l -> startswith(l, "CONECT") || startswith(l, "TITLE"), readlines(tmpfile))
end

#
# Function that returns an ATOM line in PDB format
#
# requires Printf

function _align_name(name)
    name = strip(name)
    length(name) == 1 && return " $(name)  "
    length(name) == 2 && return " $(name) "
    length(name) == 3 && return " $(name)"
    return String(name)
end
function _align_resname(resname)
    resname = strip(resname)
    length(resname) == 1 && return " $(resname) "
    length(resname) == 2 && return "$(resname)  "
    length(resname) == 3 && return "$(resname) "
    return String(resname)
end

function write_pdb_atom(atom::Atom, serial::Integer=index(atom))

    #ATOM      2  CA  GLY A   1      -1.774   6.778  32.054  1.00  0.08           C
    #COLUMNS        DATA  TYPE    FIELD        DEFINITION
    #-------------------------------------------------------------------------------------
    # 1 -  6        Record name   "ATOM  "
    # 7 - 11        Integer       serial       Atom  serial number.
    #13 - 16        Atom          name         Atom name.
    #17             Character     altLoc       Alternate location indicator.
    #18 - 20        Residue name  resName      Residue name.
    #22             Character     chainID      Chain identifier.
    #23 - 26        Integer       resSeq       Residue sequence number.
    #27             AChar         iCode        Code for insertion of residues.
    #31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    #39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    #47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    #55 - 60        Real(6.2)     occupancy    Occupancy.
    #61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    #73 - 76        String        segname      Segment identifier, left-justified (not default of PDB)
    #77 - 78        LString(2)    element      Element symbol string, right-justified.
    #79 - 80        LString(2)    charge       Charge  on the atom.

    if abs(atom.x) > 999.999 || abs(atom.y) > 999.999 || abs(atom.z) > 999.999
        println(Core.stderr,
            "Warning: coordinates of atom $(atom.index) " *
            "($(atom.name) $(atom.resname)$(atom.resnum) $(atom.chain)) " *
            "do not fit in PDB fixed format.",
        )
    end

    atom_data = (
        atom.flag == 0 ? "ATOM" : "HETATM",
        mod(serial, 1048576),
        " ",
        _align_name(atom.name),
        " ",
        _align_resname(atom.resname),
        atom.chain,
        mod(atom.resnum, 65535),
        "    ",
        atom.x,
        atom.y,
        atom.z,
        atom.occup,
        atom.beta,
        "      ",
        segname(atom),
        element(atom) === nothing ? "  " : String(element_symbol_string(atom)),
    )

    line = if serial <= 99999 && atom.resnum <= 9999 # standard integer printing
        @sprintf("%-6s%5i%1s%-4s%1s%4s%1s%4i%4s%8.3f%8.3f%8.3f%6.2f%6.2f%5s%4s%2s", atom_data...)
    elseif serial > 99999 && atom.resnum <= 9999 # Prints index in hexadecimal
        @sprintf("%-6s%5x%1s%-4s%1s%4s%1s%4i%4s%8.3f%8.3f%8.3f%6.2f%6.2f%5s%4s%2s", atom_data...)
    elseif serial <= 99999 && atom.resnum > 9999 # Prints resnum in hexadecimal code
        @sprintf("%-6s%5i%1s%-4s%1s%4s%1s%4x%4s%8.3f%8.3f%8.3f%6.2f%6.2f%5s%4s%2s", atom_data...)
    elseif serial > 99999 && atom.resnum > 9999 # Both hexadecimal
        @sprintf("%-6s%5x%1s%-4s%1s%4s%1s%4x%4s%8.3f%8.3f%8.3f%6.2f%6.2f%5s%4s%2s", atom_data...)
    end

    return line
end

@testitem "write_pdb_atom" setup = [PublicNames] begin
    using PDBToolsBase
    pdb = read_pdb(PDBToolsBase.SMALLPDB)
    @test PDBToolsBase.write_pdb_atom(pdb[1]) == "ATOM      1  N   ALA A   1      -9.229 -14.861  -5.481  0.00  0.00      PROT N"
    @test PDBToolsBase.write_pdb_atom(pdb[2]) == "ATOM      2 1HT1 ALA A   1     -10.048 -15.427  -5.569  0.00  0.00      PROT H"
    pdb[1].index = 1000000
    @test PDBToolsBase.write_pdb_atom(pdb[1]) == "ATOM  f4240  N   ALA A   1      -9.229 -14.861  -5.481  0.00  0.00      PROT N"
    pdb[1].index = 1
    pdb[1].resnum = 1000000
    @test PDBToolsBase.write_pdb_atom(pdb[1]) == "ATOM      1  N   ALA A424f      -9.229 -14.861  -5.481  0.00  0.00      PROT N"
    pdb[1].index = 1000000
    @test PDBToolsBase.write_pdb_atom(pdb[1]) == "ATOM  f4240  N   ALA A424f      -9.229 -14.861  -5.481  0.00  0.00      PROT N"
    pdb[1].resname = "ALAX"
    @test PDBToolsBase.write_pdb_atom(pdb[1]) == "ATOM  f4240  N   ALAXA424f      -9.229 -14.861  -5.481  0.00  0.00      PROT N"
    pdb[1].flag = 1
    @test PDBToolsBase.write_pdb_atom(pdb[1]) == "HETATMf4240  N   ALAXA424f      -9.229 -14.861  -5.481  0.00  0.00      PROT N"
end
