# PDBToolsBase.jl

!!! warning "Use PDBTools.jl"
    PDBToolsBase.jl is an infrastructure package. **To read, write or analyze molecular 
    structures, use [PDBTools.jl](https://m3g.github.io/PDBTools.jl/stable/)**, which 
    re-exports the public names defined here and contains the user documentation.

PDBToolsBase.jl contains the core data structures and input/output functions of PDBTools.jl,
with a small set of lightweight dependencies. It is meant to be a dependency of packages that 
need these structures but not the analysis tools of PDBTools.jl, in particular packages that 
are compiled into executables with `juliac --trim`. It also defines the function names shared
by the packages of the ecosystem (PDBTools.jl, MolSimToolkit.jl, ComplexMixtures.jl, 
Packmol.jl), previously defined in MolSimToolkitShared.jl.

```@docs
PDBToolsBase
```

The [Public API](@ref) page defines what is covered by the versioning of the package. The 
changes of each version are listed in the 
[changelog](https://github.com/m3g/PDBToolsBase.jl/blob/main/CHANGELOG.md).
