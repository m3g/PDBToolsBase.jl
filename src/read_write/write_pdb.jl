"""
    write_pdb(filename::AbstractString, atoms::AbstractVector{<:Atom}, [selection]; header=:auto, footer=:auto, append=false)

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

"""
function write_pdb(
    filename::AbstractString, 
    atoms::AbstractVector{<:Atom},
    selection_function::F=all;
    append=false,
    header=:auto,
    footer=:auto,
) where {F<:Function}
    open(expanduser(filename), append ? "a" : "w") do io
        _write_pdb_header(io, header, length(atoms))
        _write_pdb_cryst1(io, _unitcell(atoms))
        for atom in atoms
            if selection_function(atom)
                println(io, write_pdb_atom(atom))
            end
        end
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

function write_pdb_atom(atom::Atom)

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
        mod(atom.index, 1048576),
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

    line = if atom.index <= 99999 && atom.resnum <= 9999 # standard integer printing
        @sprintf("%-6s%5i%1s%-4s%1s%4s%1s%4i%4s%8.3f%8.3f%8.3f%6.2f%6.2f%5s%4s%2s", atom_data...)
    elseif atom.index > 99999 && atom.resnum <= 9999 # Prints index in hexadecimal
        @sprintf("%-6s%5x%1s%-4s%1s%4s%1s%4i%4s%8.3f%8.3f%8.3f%6.2f%6.2f%5s%4s%2s", atom_data...)
    elseif atom.index <= 99999 && atom.resnum > 9999 # Prints resnum in hexadecimal code
        @sprintf("%-6s%5i%1s%-4s%1s%4s%1s%4x%4s%8.3f%8.3f%8.3f%6.2f%6.2f%5s%4s%2s", atom_data...)
    elseif atom.index > 99999 && atom.resnum > 9999 # Both hexadecimal
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
