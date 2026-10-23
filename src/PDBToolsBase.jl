"""
    PDBToolsBase

Core data structures and input/output functions of PDBTools.jl. This is an infrastructure
package: end users should use PDBTools.jl, which re-exports everything defined here and
contains the user documentation.
"""
module PDBToolsBase

import Dates
import PrecompileTools
import Printf
using Compat: @compat
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

export atomic_number, atomic_symbol, atomic_mass, position, set_position!

export read_pdb, write_pdb, getseq, oneletter, threeletter, residuename
export read_mmcif, write_mmcif
export Atom, printatom, index, index_pdb, name, beta, occup, charge, pdb_element
export Structure, UnitCell
@compat public MetaData
export add_custom_field
export Residue, eachresidue, resname, residue, resnum, model, segname, zeta, zeta_check
export get_atoms
export Chain, eachchain, chain
export residue_ticks
export Segment, eachsegment
export Model, eachmodel
export positions, maxmin, distance, closest, dihedral, Ramachandran, set_phi!, set_psi!, extended_chain!, extended_chain
export residue_residue_distance
export element, mass, element_name, element_symbol, element_symbol_string, element_vdw_radius
export formula, stoichiometry
export Sequence

# Tools
export center_of_mass
export moveto!
export read_unitcell, lattice_to_matrix, matrix_to_lattice
export save, load

# Custom residue and element definitions
export add_protein_residue!, remove_custom_protein_residues!
export add_nucleoside_residue!, remove_custom_nucleoside_residues!
export isnucleoside, ispurine, ispyrimidine
export add_element!, remove_custom_elements!
@compat public Element, elements, ProteinResidue, protein_residues, NucleosideResidue, nucleoside_residues

# Shared function names (formerly in MolSimToolkitShared), to be extended by other packages
@compat public distances
@compat public coordination_number, bulk_coordination
@compat public wrap, wrap_to_first
@compat public align, align!, rmsd, alignment_movements, apply_alignment_transformation!
@compat public dihedrals

# Test files
const src_dir = @__DIR__
const test_dir = joinpath(src_dir, "../test")
const TESTPDB = joinpath(src_dir, "../test/structure.pdb")
const SMALLPDB = joinpath(src_dir, "../test/small.pdb")
const TESTCIF = joinpath(src_dir, "../test/small.cif")
const CHAINSPDB = joinpath(src_dir, "../test/chains.pdb")
const BROKENCIF = joinpath(src_dir, "../test/broken.cif")
const LONG_CHAIN_STRING_CIF = joinpath(src_dir, "../test/long_chain_string.cif")
const LONG_NAME_STRING_CIF = joinpath(src_dir, "../test/long_name_string.cif")
const DIMERPDB = joinpath(src_dir, "../test/dimer.pdb")
const HETATMPDB = joinpath(src_dir, "../test/hetatm.pdb")
const CIF_2C_CHAIN = joinpath(src_dir, "../test/6CO8.cif")
const TESTPBC = joinpath(src_dir, "../test/pbc_protein.pdb")
const TESTNOPBC = joinpath(src_dir, "../test/no_pbc_protein.pdb")

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

include("./properties/formula.jl")
include("./properties/getseq.jl")

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
