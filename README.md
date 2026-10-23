# PDBToolsBase

> [!IMPORTANT]
> This is an infrastructure package. **If you want to read, write or analyze molecular
> structures, use [PDBTools.jl](https://github.com/m3g/PDBTools.jl)**, which re-exports
> the public names defined here. The user documentation, including the docstrings of the
> functions defined in this package, is in the
> [PDBTools.jl documentation](https://m3g.github.io/PDBTools.jl/stable/).

PDBToolsBase contains the core data structures and input/output functions of PDBTools.jl,
with a small set of lightweight dependencies. It is meant to be a dependency of packages
that need these structures but not the analysis tools of PDBTools, in particular packages
that are compiled into executables with `juliac --trim`. It also replaces
MolSimToolkitShared.jl, defining the function names shared by the packages of the ecosystem.

The public API of each version is defined in the
[documentation](https://m3g.github.io/PDBToolsBase.jl/stable/api/), and the changes of each
version are listed in the [changelog](CHANGELOG.md).

## For developers of dependent packages

- Docstrings are written for users of PDBTools.jl: examples use `using PDBTools`, and they
  are rendered and doctested in the PDBTools.jl documentation.
- The test suite (`@testitem`s in `src/`) uses `PDBToolsBase` directly.
- `test/juliac` contains a smoke test checking that the core API compiles with
  `juliac --trim=safe`.
- Every pull request must update `CHANGELOG.md` (checked on CI; skip with the
  `Skip-Changelog` label).
