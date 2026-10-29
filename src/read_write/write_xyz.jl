"""
    write_xyz(filename::AbstractString, atoms::AbstractVector{<:Atom}, [selection]; append=false)

Write a file in the xyz format with the atoms in `atoms` to `filename`. The `selection` argument
is a string or function that can be used to select a subset of the atoms in `atoms`. For
example, `write_xyz("test.xyz", atoms, "name CA")`.

Each atom is written with its name (the `name` field), followed by its coordinates, such that
files read with [`read_xyz`](@ref) are written back with the same first column.
The comment line of each frame is the title of the structure, if `atoms` is a `Structure`
with a title, and is empty otherwise. Atoms of different models are written as different
frames.

# Keyword arguments

- `append::Bool=false`: If `true`, the atoms will be appended to the file instead of overwriting it.

"""
function write_xyz(
    filename::AbstractString,
    atoms::AbstractVector{<:Atom},
    selection_function::F=all;
    append=false,
) where {F<:Function}
    # Number of atoms of each frame: consecutive selected atoms of the same model
    frame_natoms = Int[]
    lastmodel = Int32(0)
    for atom in atoms
        selection_function(atom) || continue
        if isempty(frame_natoms) || model(atom) != lastmodel
            push!(frame_natoms, 0)
            lastmodel = model(atom)
        end
        frame_natoms[end] += 1
    end
    title = something(_title(atoms), "")
    open(expanduser(filename), append ? "a" : "w") do io
        iframe = 0
        remaining = 0 # atoms remaining to be written in the current frame
        for atom in atoms
            selection_function(atom) || continue
            if remaining == 0
                iframe += 1
                remaining = frame_natoms[iframe]
                println(io, remaining)
                println(io, title)
            end
            println(io, write_xyz_atom(atom))
            remaining -= 1
        end
    end
    return nothing
end

function write_xyz_atom(atom::Atom)
    return @sprintf("%-4s %12.5f %12.5f %12.5f", strip(name(atom)), atom.x, atom.y, atom.z)
end

@testitem "write_xyz" setup = [PublicNames] begin
    using PDBToolsBase
    s = read_xyz(PDBToolsBase.SMALLXYZ)
    tmpfile = tempname() * ".xyz"
    write_xyz(tmpfile, s)
    lines = readlines(tmpfile)
    @test lines[1:3] == ["9", "Ethanol", "C        -0.74800     -0.01500      0.02400"]
    @test length(lines) == 11
    s2 = read_xyz(tmpfile)
    @test s2 == s
    @test s2.title == "Ethanol"
    # from a PDB file, with selection
    pdb = read_pdb(joinpath(pkgdir(PDBToolsBase), "test", "ethanol.pdb"))
    write_xyz(tmpfile, pdb, at -> element(at) != "H")
    s2 = read_xyz(tmpfile)
    @test s2.title == pdb.title
    @test name.(s2) == ["C1", "C2", "O1"]
    @test positions(s2) == positions(pdb)[1:3]
    # multiple models are written as frames, and append
    s = read_xyz(PDBToolsBase.SMALLXYZ)
    s2 = Structure(vcat(s.atoms, copy.(s.atoms)))
    for at in s2[10:18]
        at.model = 2
    end
    write_xyz(tmpfile, s2)
    s3 = read_xyz(tmpfile)
    @test length(s3) == 18
    @test length(eachmodel(s3)) == 2
    @test isnothing(s3.title)
    @test readlines(tmpfile)[1:2] == ["9", ""]
    write_xyz(tmpfile, s; append=true)
    s3 = read_xyz(tmpfile)
    @test length(eachmodel(s3)) == 3
    @test positions(s3[19:27]) == positions(s)
    # the first column is written back as read: labels and atomic numbers
    xyz = "3\ncomment\nO1        0.00000      0.00000      0.00000\n" *
          "1         1.00000      0.00000      0.00000\nXY1       0.00000      1.00000      0.00000\n"
    write_xyz(tmpfile, read_xyz(IOBuffer(xyz)))
    @test read(tmpfile, String) == xyz
end
