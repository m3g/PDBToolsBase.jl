"""
    UnitCell(matrix::AbstractMatrix)
    UnitCell(sides::AbstractVector)
    UnitCell(a, b, c, α=90, β=90, γ=90)
    UnitCell(; a, b, c, α=90, β=90, γ=90)
    UnitCell(lattice::NamedTuple)

Unit cell of a periodic system. It can be defined by:

- A 3×3 matrix, with the lattice vectors as columns.
- A vector with the lengths of the sides of an orthorhombic box.
- The lattice parameters: the lengths of the sides `a`, `b`, `c` (in Å), and the angles
  `α`, `β`, `γ` (in degrees), given as positional or keyword arguments, or as the named tuple
  returned by [`matrix_to_lattice`](@ref).

A `UnitCell` behaves as the 3×3 unit cell matrix (with the lattice vectors as columns),
and can be used wherever a unit cell matrix is accepted. The unit cell of a
[`Structure`](@ref) is set with `s.unitcell = UnitCell(...)`.

# Example

```jldoctest
julia> using PDBTools

julia> uc = UnitCell(10.0, 20.0, 30.0)
UnitCell(a=10.000, b=20.000, c=30.000, α=90.00, β=90.00, γ=90.00)

julia> UnitCell([10.0, 20.0, 30.0]) == uc
true

julia> UnitCell(a=10.0, b=20.0, c=30.0, γ=120.0)
UnitCell(a=10.000, b=20.000, c=30.000, α=90.00, β=90.00, γ=120.00)

julia> matrix_to_lattice(uc)
(a = 10.0f0, b = 20.0f0, c = 30.0f0, α = 90.0f0, β = 90.0f0, γ = 90.0f0)
```

"""
struct UnitCell <: AbstractMatrix{Float32}
    matrix::SMatrix{3,3,Float32,9}
end
UnitCell(uc::UnitCell) = uc
function UnitCell(m::AbstractMatrix)
    size(m) == (3, 3) || throw(ArgumentError("The unit cell matrix must be a 3×3 matrix, got size $(size(m))."))
    return UnitCell(SMatrix{3,3,Float32,9}(m))
end
function UnitCell(sides::AbstractVector)
    length(sides) == 3 || throw(ArgumentError("The unit cell sides must be a vector of length 3, got length $(length(sides))."))
    return UnitCell(lattice_to_matrix(sides[1], sides[2], sides[3], 90, 90, 90))
end
UnitCell(a::Real, b::Real, c::Real, α::Real=90, β::Real=90, γ::Real=90) =
    UnitCell(lattice_to_matrix(Float32(a), Float32(b), Float32(c), Float32(α), Float32(β), Float32(γ)))
UnitCell(; a::Real, b::Real, c::Real, α::Real=90, β::Real=90, γ::Real=90) = UnitCell(a, b, c, α, β, γ)
UnitCell(l::NamedTuple{(:a, :b, :c, :α, :β, :γ)}) = UnitCell(l.a, l.b, l.c, l.α, l.β, l.γ)

# Matrix interface
Base.size(::UnitCell) = (3, 3)
Base.@propagate_inbounds Base.getindex(uc::UnitCell, i::Int) = uc.matrix[i]
Base.@propagate_inbounds Base.getindex(uc::UnitCell, i::Int, j::Int) = uc.matrix[i, j]
Base.IndexStyle(::Type{UnitCell}) = IndexLinear()

# Specialized methods using the static matrix
matrix_to_lattice(uc::UnitCell) = matrix_to_lattice(uc.matrix)
wrap(x::AbstractVector, xref::AbstractVector, uc::UnitCell) = wrap(x, xref, uc.matrix)
wrap_to_first(x::AbstractVector{<:AbstractFloat}, uc::UnitCell) = wrap_to_first(x, uc.matrix)

function Base.show(io::IO, uc::UnitCell)
    l = matrix_to_lattice(uc)
    print(io, "UnitCell", @sprintf("(a=%.3f, b=%.3f, c=%.3f, α=%.2f, β=%.2f, γ=%.2f)", l.a, l.b, l.c, l.α, l.β, l.γ))
end
Base.show(io::IO, ::MIME"text/plain", uc::UnitCell) = show(io, uc)

@testitem "UnitCell" setup = [PublicNames] begin
    using PDBToolsBase
    using StaticArrays: SMatrix
    uc = UnitCell(10.0, 20.0, 30.0)
    @test uc isa AbstractMatrix{Float32}
    @test uc.matrix isa SMatrix{3,3,Float32,9}
    @test uc ≈ [10 0 0; 0 20 0; 0 0 30]
    @test size(uc) == (3, 3)
    @test UnitCell([10, 20, 30]) == uc
    @test UnitCell([10 0 0; 0 20 0; 0 0 30]) == uc
    @test UnitCell(uc) === uc
    @test UnitCell(a=10, b=20, c=30) == uc
    @test UnitCell(matrix_to_lattice(uc)) ≈ uc
    l = (a=39.85, b=54.45, c=59.84, α=72.01, β=89.44, γ=89.43)
    tric = UnitCell(; l...)
    @test tric ≈ lattice_to_matrix(l...)
    @test all(isapprox.(values(matrix_to_lattice(tric)), values(l); rtol=1e-5))
    @test UnitCell(l) == tric
    @test_throws ArgumentError UnitCell([1 0; 0 1])
    @test_throws ArgumentError UnitCell([1, 2])
    @test_throws UndefKeywordError UnitCell(a=1, b=2)
    # wrapping
    @test PDBToolsBase.wrap([11.0, 1.0, 1.0], [0.0, 0.0, 0.0], uc) ≈ [1.0, 1.0, 1.0]
    @test PDBToolsBase.wrap_to_first([11.0, 21.0, 31.0], uc) ≈ [1.0, 1.0, 1.0] atol = 1e-4
    @test UnitCell([10 0 0; 0 20 0; 0 0 30]).matrix == [10 0 0; 0 20 0; 0 0 30]
    # show
    @test sprint(show, uc) == "UnitCell(a=10.000, b=20.000, c=30.000, α=90.00, β=90.00, γ=90.00)"
    @test sprint(show, MIME"text/plain"(), uc) == sprint(show, uc)
end
