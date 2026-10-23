#
# Smoke test for `juliac --trim=safe` compatibility of the core API. From this directory:
#
#   JULIA_APPS_JULIA_CMD=$(which julia) juliac --output-exe pdbtoolsbase_smoke --project . \
#       --experimental --trim=safe --bundle build app.jl
#   ./build/bin/pdbtoolsbase_smoke ../structure.pdb output
#
# (Julia 1.13+). Printing uses `Core.stdout`, because `stdout` is not initialized in
# trimmed executables.
#
using PDBToolsBase: read_pdb, read_mmcif, write_pdb, write_mmcif, read_unitcell,
    isprotein, name, eachresidue, eachchain, center_of_mass, moveto!, Structure, UnitCell

function (@main)(ARGS)
    input, output = ARGS[1], ARGS[2]
    atoms = read_pdb(input)
    protein = read_pdb(input, isprotein)
    ca = filter(at -> name(at) == "CA", protein)
    nres = length(eachresidue(protein))
    nchains = length(eachchain(protein))
    cm = center_of_mass(protein)
    moveto!(protein)
    write_pdb(output * ".pdb", protein)
    write_mmcif(output * ".cif", protein)
    back_pdb = read_pdb(output * ".pdb")
    back_cif = read_mmcif(output * ".cif")
    # Structure with unit cell
    s = Structure(protein.atoms; unitcell=UnitCell(80.0, 80.0, 80.0))
    s.unitcell = UnitCell(a=80.0, b=80.0, c=80.0, γ=90.0)
    s_ca = filter(at -> name(at) == "CA", s)
    s_view = @view s[1:10]
    write_pdb(output * "_structure.pdb", s_ca)
    write_mmcif(output * "_structure.cif", s_view)
    uc_pdb = read_unitcell(output * "_structure.pdb")
    uc_cif = read_unitcell(output * "_structure.cif")
    println(Core.stdout, "atoms: ", length(atoms))
    println(Core.stdout, "protein atoms: ", length(protein), " CA: ", length(ca))
    println(Core.stdout, "residues: ", nres, " chains: ", nchains)
    println(Core.stdout, "center of mass: ", cm[1], " ", cm[2], " ", cm[3])
    println(Core.stdout, "read back: ", length(back_pdb), " ", length(back_cif))
    println(Core.stdout, "structure CA: ", length(s_ca), " unit cell read back: ", uc_pdb[1, 1], " ", uc_cif[3, 3])
    println(Core.stdout, "unit cell of input: ", isnothing(atoms.unitcell) ? "none" : sprint(show, atoms.unitcell))
    println(Core.stdout, "unit cell of written structure: ", sprint(show, read_pdb(output * "_structure.pdb").unitcell))
    return 0
end
