"""
    PDBToolsBase

Core data structures and input/output functions of PDBTools.jl. This is an infrastructure
package: end users should use PDBTools.jl, which re-exports the public names defined here and
contains the user documentation.
"""
module PDBToolsBase

import Dates
import PrecompileTools
import Printf
using LinearAlgebra: eigvecs, norm, cross, dot
using Printf: @printf, @sprintf
using StaticArrays: SVector, SMatrix, MMatrix
using TestItems: @testitem, @testmodule
using OrderedCollections: OrderedDict
using InlineStrings: InlineString, String1, String15

# Version
const VERSION = pkgversion(@__MODULE__)

# Default fixed String size for names, chains, segments, etc.
const StringType = String15

# The package exports nothing: the names below are public, and are imported explicitly
# by the packages that use them (PDBTools.jl re-exports them)
public atomic_number, atomic_symbol, atomic_mass, position, set_position!

public read_pdb, write_pdb, oneletter, threeletter, residuename
public read_mmcif, write_mmcif
public read_xyz, write_xyz
public Atom, printatom, index, index_pdb, name, beta, occup, charge, pdb_element
public Structure, UnitCell, MetaData
public add_custom_field
public Residue, eachresidue, resname, residue, resnum, model, segname
public get_atoms
public Chain, eachchain, chain
public Segment, eachsegment
public Model, eachmodel
public positions, maxmin, distance, closest, dihedral
public residue_residue_distance
public element, mass, element_name, element_symbol, element_symbol_string, element_vdw_radius
public Sequence

# Classification of atoms and residues
public isprotein, isnucleoside, isbackbone, issidechain, iswater
public isacidic, isaliphatic, isaromatic, isbasic, ischarged
public ishydrophobic, isneutral, isnonpolar, ispolar
public ispurine, ispyrimidine

# Tools
public center_of_mass
public moveto!
public read_unitcell, lattice_to_matrix, matrix_to_lattice
public save, load

# Custom residue and element definitions
public add_protein_residue!, remove_custom_protein_residues!
public add_nucleoside_residue!, remove_custom_nucleoside_residues!
public add_element!, remove_custom_elements!
public Element, elements, ProteinResidue, protein_residues, NucleosideResidue, nucleoside_residues

# Shared function names (formerly in MolSimToolkitShared), to be extended by other packages
public distances
public coordination_number
public wrap, wrap_to_first
public align, align!, rmsd, alignment_movements, apply_alignment_transformation!
public dihedrals

# Test files read by the precompilation workload, which are part of the package
const src_dir = @__DIR__
const SMALLPDB = joinpath(src_dir, "../test/small.pdb")
const TESTCIF = joinpath(src_dir, "../test/small.cif")
const SMALLXYZ = joinpath(src_dir, "../test/small.xyz")

# test_file(name): path of the test file `name`. The test files are a lazy artifact (see
# Artifacts.toml), downloaded on the first call. The method is defined in the TestFilesExt
# extension, loaded with `using LazyArtifacts`.
function test_file end

# Function that returns true for all atoms: the default selection
all(atoms) = true

#
# Shared function names and generic utilities (formerly in MolSimToolkitShared)
#
include("./shared/shared_names.jl")
include("./shared/wrap.jl")
include("./shared/structural_alignment.jl")
include("./shared/dihedral.jl")

# Basic chemistry
include("./properties/elements.jl")
include("./properties/protein_residues.jl")
include("./properties/nucleoside_residues.jl")
include("./properties/residue_properties.jl")

#
# Data structures and their iterators
#
include("./struct_and_iterators/Atom.jl")
include("./struct_and_iterators/UnitCell.jl")
include("./struct_and_iterators/Structure.jl")
include("./struct_and_iterators/iterators_base.jl")
include("./struct_and_iterators/Residue.jl")
include("./struct_and_iterators/Chain.jl")
include("./struct_and_iterators/Segment.jl")
include("./struct_and_iterators/Model.jl")

#
# Input and output functions
#
include("./read_write/parsers.jl")
include("./read_write/read_pdb.jl")
include("./read_write/write_pdb.jl")
include("./read_write/read_mmcif.jl")
include("./read_write/write_mmcif.jl")
include("./read_write/read_xyz.jl")
include("./read_write/write_xyz.jl")


# Coordinates-dependent functions
include("./coordinates/positions.jl")
include("./coordinates/distance.jl")
include("./coordinates/maxmin.jl")
include("./coordinates/dihedrals.jl")

# Miscellaneous tools
include("./tools/center_of_mass.jl")
include("./tools/read_unitcell.jl")

# Precompilation
include("./precompile.jl")

end
