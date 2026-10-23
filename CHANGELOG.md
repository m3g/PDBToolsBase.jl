PDBToolsBase.jl Changelog
===========================

[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-experimental]: https://img.shields.io/badge/Experimental-yellow.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-fix]: https://img.shields.io/badge/Fix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg

Version 1.0.0-DEV
--------------
- ![INFO][badge-info] First version. PDBToolsBase.jl contains the core data structures and input/output functions of PDBTools.jl (from PDBTools 3.41), and the shared function names and utilities of MolSimToolkitShared.jl (from MolSimToolkitShared 1.6), which it replaces. The public API is defined in the documentation.
- ![INFO][badge-info] The package exports nothing: its API is formed by the names declared `public`, which are used qualified (`PDBToolsBase.read_pdb`) or imported explicitly (`using PDBToolsBase: read_pdb`). Requires Julia ≥ 1.11.
- ![INFO][badge-info] Analysis and modelling functions remain defined in PDBTools.jl and are not part of PDBToolsBase.jl: `residue_ticks`, `formula`, `stoichiometry`, `getseq`, `zeta`, `zeta_check`, `Ramachandran`, `set_phi!`, `set_psi!`, `extended_chain` and `extended_chain!`. The shared name `bulk_coordination` of MolSimToolkitShared.jl is not carried over, as only MolSimToolkit.jl defines methods for it.
- ![FEATURE][badge-feature] `Structure` type, which behaves as a vector of atoms and carries metadata (currently, the unit cell). `read_pdb` and `read_mmcif` return a `Structure`, with the unit cell read from the `CRYST1` (PDB) or `_cell` (mmCIF) fields. `write_pdb` and `write_mmcif` write the unit cell of `Structure`s.
- ![FEATURE][badge-feature] `UnitCell` type, constructed from a matrix, a vector of box sides, or lattice parameters, and set with `s.unitcell = UnitCell(...)`.
- ![FEATURE][badge-feature] Functions that support periodic boundary conditions accept `unitcell=true` to use the unit cell of the `Structure` (`false`, the default, for no periodicity; or a unit cell).
- ![FEATURE][badge-feature] The core API (reading, writing and manipulating structures) compiles into executables with `juliac --trim=safe` (Julia ≥ 1.13).
- ![ENHANCEMENT][badge-enhancement] Type-stable and faster reading and writing of PDB and mmCIF files (`read_pdb` with 4x fewer allocations; `write_mmcif` 1.6x faster). The mmCIF reader no longer creates types from the column layout of each file.
- ![ENHANCEMENT][badge-enhancement] `lattice_to_matrix` accepts integer arguments, and gives exact zeros for right angles.
- ![FIX][badge-fix] The warning for coordinates that do not fit the PDB format checked the `y` coordinate twice, and not the `z` coordinate.
- ![FIX][badge-fix] The `write_mmcif` docstring was attached to an internal constant.
