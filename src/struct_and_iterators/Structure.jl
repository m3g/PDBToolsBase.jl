"""
    MetaData

Metadata of a [`Structure`](@ref). The fields of `MetaData` are accessible as properties
of the structure (e.g. `s.unitcell`), and can be set in the `Structure` constructor.

# Fields

- `unitcell::Union{Nothing,UnitCell}`: The unit cell of the system (see [`UnitCell`](@ref)), or 
  `nothing` if the structure is not periodic. When setting it (`s.unitcell = ...`), values that are 
  not a `UnitCell` are converted with the `UnitCell` constructor.

!!! note
    All fields of `MetaData` accept `nothing`, which means that the property is not defined.
    When structures with different metadata are concatenated, the conflicting properties are
    set to `nothing`.

"""
mutable struct MetaData
    unitcell::Union{Nothing,UnitCell}
end
MetaData(; unitcell=nothing) = MetaData(_unitcell_value(unitcell))

_unitcell_value(::Nothing) = nothing
_unitcell_value(unitcell) = UnitCell(unitcell)

function Base.setproperty!(m::MetaData, field::Symbol, value)
    field === :unitcell && return setfield!(m, :unitcell, _unitcell_value(value))
    return setfield!(m, field, value)
end
Base.copy(m::MetaData) = MetaData((getfield(m, f) for f in fieldnames(MetaData))...)

# Metadata of concatenated structures: properties that differ are set to `nothing`
function _merge_metadata(m1::MetaData, ms::MetaData...)
    values = ntuple(fieldcount(MetaData)) do i
        v = getfield(m1, i)
        Base.all(m -> isequal(getfield(m, i), v), ms) ? v : nothing
    end
    return MetaData(values...)
end

"""
    Structure(atoms::AbstractVector{<:Atom}=Atom{Nothing}[]; unitcell=nothing)
    Structure(atoms::AbstractVector{<:Atom}, metadata::MetaData)

A `Structure` is the main data structure in PDBTools.jl representing a molecular structure.
It behaves as a vector of `Atom` objects, and carries metadata (see `PDBTools.MetaData`),
such as the unit cell of the system, which is accessible as a property: `s.unitcell`.
The underlying vector of atoms is `s.atoms`.

Indexing with ranges, `filter`, `copy` and other functions that create new vectors of atoms
return new `Structure`s, with copies of the metadata, while `view` returns a `Structure`
sharing the metadata with the original one.

# Example

```jldoctest
julia> using PDBTools

julia> s = Structure(read_pdb(PDBTools.test_file("small.pdb")); unitcell=UnitCell(10.0, 10.0, 10.0));

julia> s.unitcell
UnitCell(a=10.000, b=10.000, c=10.000, α=90.00, β=90.00, γ=90.00)

julia> s.unitcell = UnitCell(a=10.0, b=20.0, c=30.0, γ=120.0)
UnitCell(a=10.000, b=20.000, c=30.000, α=90.00, β=90.00, γ=120.00)

julia> length(s[1:10])
10

julia> s[1:10].unitcell == s.unitcell
true
```

"""
struct Structure{T<:Atom,V<:AbstractVector{T}} <: AbstractVector{T}
    atoms::V
    metadata::MetaData
end
Structure(atoms::AbstractVector{<:Atom}=Atom{Nothing}[]; kargs...) = Structure(atoms, MetaData(; kargs...))

# Unit cell of a vector of atoms, if available
_unitcell(::AbstractVector{<:Atom}) = nothing
_unitcell(s::Structure) = s.unitcell

#
# Interpretation of the `unitcell` keyword argument of functions that support periodic
# boundary conditions: `false` (or `nothing`): no periodicity; `true`: use the unit cell of
# the structure; a matrix or vector: the unit cell to be used. Unit cells of type `UnitCell` 
# are returned as static matrices.
#
_resolve_unitcell(::AbstractVector{<:Atom}, ::Nothing) = nothing
_resolve_unitcell(::AbstractVector{<:Atom}, unitcell::AbstractVecOrMat) = unitcell
_resolve_unitcell(::AbstractVector{<:Atom}, unitcell::UnitCell) = unitcell.matrix
function _resolve_unitcell(atoms::AbstractVector{<:Atom}, unitcell::Bool)
    unitcell || return nothing
    uc = _unitcell(atoms)
    if isnothing(uc)
        throw(ArgumentError("""\n
            unitcell=true requires a Structure with a unit cell, but the unit cell is not defined.
            Set it with `s.unitcell = ...`, or provide the unit cell matrix with `unitcell=...`.

        """))
    end
    return uc.matrix
end

#
# Properties: the atoms, the metadata, and the fields of the metadata
#
function Base.getproperty(s::Structure, field::Symbol)
    field === :atoms && return getfield(s, :atoms)
    field === :metadata && return getfield(s, :metadata)
    return getproperty(getfield(s, :metadata), field)
end
function Base.setproperty!(s::Structure, field::Symbol, value)
    if field === :atoms || field === :metadata
        throw(ArgumentError("Cannot replace the $field of a Structure."))
    end
    return setproperty!(getfield(s, :metadata), field, value)
end
Base.propertynames(::Structure) = (:atoms, :metadata, fieldnames(MetaData)...)

#
# Array interface
#
Base.size(s::Structure) = size(getfield(s, :atoms))
Base.IndexStyle(::Type{<:Structure}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(s::Structure, i::Int) = getfield(s, :atoms)[i]
Base.@propagate_inbounds Base.setindex!(s::Structure, at::Atom, i::Int) = setindex!(getfield(s, :atoms), at, i)

# `similar` is used by `getindex` with ranges, `sort`, `map`, etc. Vectors of atoms are Structures
# with copies of the metadata, other vectors are not.
Base.similar(s::Structure, ::Type{T}, dims::Dims{1}) where {T<:Atom} =
    Structure(similar(getfield(s, :atoms), T, dims), copy(getfield(s, :metadata)))
Base.similar(s::Structure, ::Type{T}, dims::Dims) where {T} = similar(getfield(s, :atoms), T, dims)
Base.view(s::Structure, inds::Union{Colon,AbstractVector{<:Integer}}) =
    Structure(view(getfield(s, :atoms), inds), getfield(s, :metadata))
Base.copy(s::Structure) = Structure(copy(getfield(s, :atoms)), copy(getfield(s, :metadata)))
Base.filter(f, s::Structure) = Structure(filter(f, getfield(s, :atoms)), copy(getfield(s, :metadata)))
Base.vcat(s::Structure, ss::Structure...) = Structure(
    vcat(getfield(s, :atoms), (getfield(x, :atoms) for x in ss)...),
    _merge_metadata(getfield(s, :metadata), (getfield(x, :metadata) for x in ss)...),
)

# Mutating functions act on the underlying vector of atoms
Base.push!(s::Structure, at::Atom) = (push!(getfield(s, :atoms), at); s)
Base.pushfirst!(s::Structure, at::Atom) = (pushfirst!(getfield(s, :atoms), at); s)
Base.insert!(s::Structure, i::Integer, at::Atom) = (insert!(getfield(s, :atoms), i, at); s)
Base.append!(s::Structure, ats) = (append!(getfield(s, :atoms), ats); s)
Base.prepend!(s::Structure, ats) = (prepend!(getfield(s, :atoms), ats); s)
Base.resize!(s::Structure, n::Integer) = (resize!(getfield(s, :atoms), n); s)
Base.sizehint!(s::Structure, n::Integer) = (sizehint!(getfield(s, :atoms), n); s)
Base.deleteat!(s::Structure, inds) = (deleteat!(getfield(s, :atoms), inds); s)
Base.keepat!(s::Structure, inds) = (keepat!(getfield(s, :atoms), inds); s)
Base.filter!(f, s::Structure) = (filter!(f, getfield(s, :atoms)); s)
Base.empty!(s::Structure) = (empty!(getfield(s, :atoms)); s)
Base.pop!(s::Structure) = pop!(getfield(s, :atoms))
Base.popfirst!(s::Structure) = popfirst!(getfield(s, :atoms))

function Base.show(io::IO, s::Structure{T}) where {T}
    lines, cols = _displaysize(io)
    compact = get(io, :compact, false)::Bool
    ntitlelines = 5
    type = get(io, :type, true)::Bool
    if type && !compact && cols >= 115 && lines > ntitlelines - 1
        cell = isnothing(s.unitcell) ? "" : " - with unit cell"
        println(io, "   Structure{$T} with $(length(s)) atoms$cell")
    end
    show(IOContext(io, :type => false, :ntitlelines => ntitlelines), s.atoms)
end

@testitem "Structure" setup = [PublicNames] begin
    using PDBToolsBase
    using StaticArrays: SMatrix
    s = Structure(read_pdb(PDBToolsBase.test_file("structure.pdb")))
    ats = s.atoms
    @test s[1] == ats[1]
    @test s[end] == ats[end]
    @test s[5:10] == ats[5:10]
    @test size(s) == size(ats)
    @test length(s) == length(ats)
    s[2] = Atom(name="CA")
    @test s[2] == Atom(name="CA")
    @test propertynames(s) == (:atoms, :metadata, :unitcell)
    @test isnothing(s.unitcell)
    s.unitcell = [ 1 0 0 ; 0 1 0 ; 0 0 1]
    @test s.unitcell == [ 1 0 0 ; 0 1 0 ; 0 0 1]
    @test s.unitcell isa UnitCell
    s.unitcell = [10, 20, 30]
    @test s.unitcell == UnitCell(10, 20, 30)
    s.unitcell = UnitCell(a=10, b=20, c=30, γ=120)
    @test matrix_to_lattice(s.unitcell).γ ≈ 120
    s.unitcell = nothing
    @test isnothing(s.unitcell)
    s = Structure(ats; unitcell=[1 0 0; 0 1 0; 0 0 1])
    @test s.unitcell == [ 1 0 0 ; 0 1 0 ; 0 0 1]
    @test_throws "Cannot replace" s.metadata = PDBToolsBase.MetaData()
    @test_throws "Cannot replace" s.atoms = ats
    @test_throws ArgumentError s.unitcell = [1 0; 0 1]
    @test_throws ArgumentError Structure(ats; unitcell=[1 0; 0 1])
    @test_throws MethodError Structure(ats; abc=1)
end

@testitem "Structure - vector interface" setup = [PublicNames] begin
    using PDBToolsBase
    using StaticArrays: SMatrix
    uc = [10 0 0; 0 10 0; 0 0 10]
    s = Structure(read_pdb(PDBToolsBase.SMALLPDB); unitcell=uc)
    n = length(s)
    # slices are copies, with copies of the metadata
    s2 = s[1:10]
    @test s2 isa Structure
    @test s2.unitcell == uc
    @test s2.atoms !== s.atoms
    s2.unitcell = nothing
    @test s.unitcell == uc
    @test s[[1, 3]] == [s[1], s[3]]
    @test s[[true; falses(n - 1)]] == [s[1]]
    # views share the atoms and the metadata
    v = @view s[2:5]
    @test v isa Structure
    @test v[1] === s[2]
    v.unitcell = 2 .* uc
    @test s.unitcell == 2 .* uc
    s.unitcell = uc
    # functions that return vectors of atoms return Structures
    @test copy(s) isa Structure
    @test copy(s).unitcell == uc
    @test copy(s).atoms !== s.atoms
    sel = filter(isprotein, s)
    @test sel isa Structure
    @test sel.unitcell == uc
    @test length(sel) == count(isprotein, s)
    @test reverse(s) isa Structure
    @test sort(s; by=at -> -index(at))[1] === s[end]
    @test map(identity, s) isa Structure
    @test map(position, s) isa Vector
    @test position.(s) isa Vector
    @test positions(s) == positions(s.atoms)
    # concatenation keeps the metadata only if it is the same
    @test vcat(s, s) isa Structure
    @test length(vcat(s, s)) == 2n
    @test vcat(s, s).unitcell == uc
    @test isnothing(vcat(s, Structure(s.atoms)).unitcell)
    # mutating functions
    s3 = copy(s)
    push!(s3, Atom(name="X1"))
    @test length(s3) == n + 1 && name(s3[end]) == "X1"
    pushfirst!(s3, Atom(name="X0"))
    @test name(s3[1]) == "X0"
    insert!(s3, 2, Atom(name="X2"))
    @test name(s3[2]) == "X2"
    @test name(pop!(s3)) == "X1"
    @test name(popfirst!(s3)) == "X0"
    deleteat!(s3, 1)
    @test s3 == s
    append!(s3, s)
    @test length(s3) == 2n
    resize!(s3, n)
    @test s3 == s
    prepend!(s3, [Atom(name="X3")])
    @test name(s3[1]) == "X3"
    keepat!(s3, 2:n+1)
    @test s3 == s
    filter!(isprotein, s3)
    @test length(s3) == count(isprotein, s)
    @test empty!(s3) isa Structure
    @test isempty(s3)
    @test s3.unitcell == uc
    # iterators and functions on vectors of atoms
    @test length(eachresidue(s)) == length(eachresidue(s.atoms))
    @test center_of_mass(s) ≈ center_of_mass(s.atoms)
    @test PDBToolsBase._unitcell(s) == uc
    @test isnothing(PDBToolsBase._unitcell(s.atoms))
    @test s isa Structure{Atom{Nothing}}
    # interpretation of the unitcell keyword
    @test PDBToolsBase._resolve_unitcell(s, true) == uc
    @test PDBToolsBase._resolve_unitcell(s, true) isa SMatrix{3,3,Float32,9}
    @test PDBToolsBase._resolve_unitcell(s, UnitCell(uc)) isa SMatrix{3,3,Float32,9}
    @test isnothing(PDBToolsBase._resolve_unitcell(s, false))
    @test isnothing(PDBToolsBase._resolve_unitcell(s, nothing))
    @test PDBToolsBase._resolve_unitcell(s, [1.0, 2.0, 3.0]) == [1.0, 2.0, 3.0]
    @test PDBToolsBase._resolve_unitcell(s.atoms, 2 .* uc) == 2 .* uc
    @test_throws "unitcell=true requires" PDBToolsBase._resolve_unitcell(s.atoms, true)
    @test_throws "unitcell=true requires" PDBToolsBase._resolve_unitcell(Structure(s.atoms), true)
end

@testitem "Structure - show" setup = [PublicNames] begin
    using PDBToolsBase
    using ShowMethodTesting
    ENV["LINES"] = 10
    ENV["COLUMNS"] = 120
    at = Atom(name="CA", segname="X")
    s = Structure([at]; unitcell=[1 0 0; 0 1 0; 0 0 1])
    @test parse_show(s; repl=Dict("PDBToolsBase." => "")) ≈ """
       Structure{Atom{Nothing}} with 1 atoms - with unit cell
       index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       0   CA     XXX     X        0        0    0.000    0.000    0.000  0.00  0.00     0       X         0
    """
    @test parse_show(Structure([at for _ in 1:50]); repl=Dict(r"^((?:[^\n]*\n){4}).*"s => s"\1", r"PDBToolsBase." => "")) ≈ """
       Structure{Atom{Nothing}} with 50 atoms
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       0   CA     XXX     X        0        0    0.000    0.000    0.000  0.00  0.00     0       X         0
    """
    @test parse_show(Dict(1 => Structure([at, at, at])); repl=Dict("PDBToolsBase." => "")) ≈ """
    Dict{Int64, Structure{Atom{Nothing}, Vector{Atom{Nothing}}}} with 1 entry:
    1 => [ Atom(0CA-XXX0X), Atom(0CA-XXX0X), Atom(0CA-XXX0X) ]
    """
    @test parse_show([Structure([at, at]), Structure([at, at, at])]; repl=Dict("PDBToolsBase." => "")) ≈ """
    2-element Vector{Structure{Atom{Nothing}, Vector{Atom{Nothing}}}}[
    [ Atom(0CA-XXX0X), Atom(0CA-XXX0X) ]
    [ Atom(0CA-XXX0X), Atom(0CA-XXX0X), Atom(0CA-XXX0X) ]
    ]
    """
    @test parse_show(Atom{Nothing}[]; repl=Dict("PDBToolsBase." => "")) ≈ """
       Vector{Atom{Nothing}} with 0 atoms with fields:
    index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
    """
    @test parse_show(Structure(); repl=Dict("PDBToolsBase." => "")) ≈ """
       Structure{Atom{Nothing}} with 0 atoms
    index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
    """
    @test parse_show([ Structure() ]; repl=Dict("PDBToolsBase." => "")) ≈ """
    1-element Vector{Structure{Atom{Nothing}, Vector{Atom{Nothing}}}}[
        [  ]
    ]
    """
    @test parse_show(Structure([ at for _ in 1:20 ]); mime=nothing, context = :compact => true) ≈ """
    [ Atom( 0 CA-XXX 0 X), Atom( 0 CA-XXX 0 X), Atom( 0 CA-XXX 0 X), Atom( 0 CA-XXX 0 X)…
    """
end
