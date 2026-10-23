# Public API

PDBToolsBase.jl follows [semantic versioning](https://semver.org/). This page defines the 
public API of the version of the package it belongs to: what is listed here is only changed
in a breaking way in a new major version. Anything not listed here is internal, and may 
change in any release.

!!! note
    The docstrings are written for users of PDBTools.jl, which re-exports this package. 
    Thus, the examples use `using PDBTools`, and refer to `PDBTools.name` for names that are
    public but not exported.

## What is public

1. **Public names**: the names exported by `PDBToolsBase`, and the names declared `public` 
   (`Base.ispublic(PDBToolsBase, name) == true`, in Julia ≥ 1.11). All of them are 
   documented in this page (this is checked when the documentation is built), with the 
   methods and keyword arguments described in their docstrings.

2. **Documented fields and properties**:
    - The fields of [`Atom`](@ref): `index`, `index_pdb`, `name`, `resname`, `chain`, 
      `resnum`, `residue`, `x`, `y`, `z`, `beta`, `occup`, `model`, `segname`, 
      `pdb_element`, `charge`, and `custom` (see [`add_custom_field`](@ref)).
    - The properties of [`Structure`](@ref): `atoms`, `metadata`, and the fields of 
      [`PDBToolsBase.MetaData`](@ref) (currently, `unitcell`).
    - The `matrix` field of [`UnitCell`](@ref).
    - The fields of [`PDBToolsBase.Element`](@ref), [`PDBToolsBase.ProteinResidue`](@ref), and
      [`PDBToolsBase.NucleosideResidue`](@ref).

3. **Documented behavior**:
    - [`read_pdb`](@ref) and [`read_mmcif`](@ref) return a `Structure`, with the unit cell
      read from the file (`CRYST1` or `_cell` fields), or `nothing` if the file does not 
      define a unit cell or defines the 1Å unit cube (which indicates a non-crystallographic
      structure).
    - Selections are functions that receive an `Atom` and return `true` or `false` (e.g. 
      `isprotein`, or `at -> name(at) == "CA"`). Selections given as strings 
      (`"protein and name CA"`) are provided by PDBTools.jl.
    - Functions that support periodic boundary conditions accept the `unitcell` keyword
      argument as: `false` (no periodicity), `true` (use the unit cell of the `Structure`), 
      or a unit cell (a `UnitCell`, a 3×3 matrix, or a vector of box sides).
    - A `Structure` behaves as a vector of atoms (`AbstractVector{<:Atom}`). Indexing with
      ranges, `filter`, `copy`, and other functions that create new vectors of atoms 
      return new `Structure`s with copies of the metadata, while `view` returns a 
      `Structure` that shares the atoms and the metadata with the original one. When 
      structures are concatenated, metadata that differs between them is set to `nothing`.
    - The files written by [`write_pdb`](@ref) and [`write_mmcif`](@ref) follow the PDB and
      PDBx/mmCIF formats, including the unit cell of `Structure`s.
    - The core API (reading, writing, and manipulating structures) can be compiled into 
      executables with `juliac --trim=safe` (Julia ≥ 1.13).

## What is not public

- Names that are not exported or declared `public`, even if accessible as 
  `PDBToolsBase.name` (for example, the test files `PDBToolsBase.TESTPDB`, or internal 
  functions with names starting with `_`).
- Methods and keyword arguments that are not described in the docstrings.
- The printed output (`show`) of the types, and the exact formatting (column widths, 
  number of digits, header lines) of the files written, as long as they remain valid PDB
  or mmCIF files.
- The `flag` field of `Atom`.
- The values of the properties of elements and residues, which may be corrected in any
  release.

Packages of the ecosystem developed together with PDBToolsBase.jl (in particular PDBTools.jl)
may use internal names, and declare compatibility bounds accordingly. Other packages should
use only the public API.

## Data structures

```@docs
Atom
Structure
PDBToolsBase.MetaData
UnitCell
Residue
Chain
Segment
Model
Sequence
```

## Reading and writing

```@docs
read_pdb
read_mmcif
write_pdb
write_mmcif
```

## Atom properties

```@docs
index
index_pdb
name
resname
chain
resnum
residue
segname
model
beta
occup
charge
pdb_element
PDBToolsBase.position(::Atom)
set_position!
printatom
add_custom_field
get_atoms
```

## Iterators

```@docs
eachresidue
eachchain
eachsegment
eachmodel
residue_ticks
```

## Elements and residues

```@docs
element
element_name
element_symbol
element_symbol_string
element_vdw_radius
atomic_number
atomic_symbol
atomic_mass
mass
formula
stoichiometry
getseq
oneletter
threeletter
residuename
zeta
zeta_check
```

## Classification of atoms and residues

```@docs
isprotein
isnucleoside
isbackbone
issidechain
iswater
isacidic
isaliphatic
isaromatic
isbasic
ischarged
isneutral
ishydrophobic
ispolar
isnonpolar
ispurine
ispyrimidine
```

## Custom elements and residues

```@docs
add_element!
remove_custom_elements!
add_protein_residue!
remove_custom_protein_residues!
add_nucleoside_residue!
remove_custom_nucleoside_residues!
PDBToolsBase.Element
PDBToolsBase.elements
PDBToolsBase.ProteinResidue
PDBToolsBase.protein_residues
PDBToolsBase.NucleosideResidue
PDBToolsBase.nucleoside_residues
```

## Coordinates

```@docs
positions
distance
closest
residue_residue_distance
maxmin
center_of_mass
moveto!
```

## Dihedral angles

```@docs
dihedral
PDBToolsBase.dihedrals
Ramachandran
set_phi!
set_psi!
extended_chain
extended_chain!
```

## Periodic systems

```@docs
lattice_to_matrix
matrix_to_lattice
read_unitcell
PDBToolsBase.wrap
PDBToolsBase.wrap_to_first
```

## Structural alignment

```@docs
PDBToolsBase.align
PDBToolsBase.align!
PDBToolsBase.rmsd
PDBToolsBase.alignment_movements
PDBToolsBase.apply_alignment_transformation!
```

## Shared function names

These functions have methods defined in other packages of the ecosystem, which extend the 
functions defined here.

```@docs
PDBToolsBase.distances
PDBToolsBase.bulk_coordination
PDBToolsBase.coordination_number
save
load
```
