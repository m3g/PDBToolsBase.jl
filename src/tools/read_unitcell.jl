"""
    read_unitcell(file::AbstractString)

Reads the lattice parameters of unitcell from a PDB (`CRYST1` field) or mmCIF file, and converts it to a unitcell matrix
with the `lattice_to_matrix` function. The unitcell matrix contains the lattice vectors as the *columns* of the matrix. 

# Example

```julia-repl
julia> using PDBTools

julia> m = read_unitcell(PDBTools.TESTPBC)
3×3 StaticArraysCore.SMatrix{3, 3, Float32, 9} with indices SOneTo(3)×SOneTo(3):
 85.0   0.0   0.0
  0.0  85.0   0.0
  0.0   0.0  85.0

julia> matrix_to_lattice(m)
(a = 85.0f0, b = 85.0f0, c = 85.0f0, α = 90.0f0, β = 90.0f0, γ = 90.0f0)
```

"""
function read_unitcell(file::AbstractString)
    lattice = open(_read_lattice_parameters, file)
    if !_is_complete(lattice)
        throw(ArgumentError("Could not find (complete?) unit cell information in file."))
    end
    return lattice_to_matrix(lattice...)
end

# Lattice parameters (a, b, c, α, β, γ): the parameters not available are `NaN`
# (a concrete type is required for `juliac --trim`)
const _Lattice = NTuple{6,Float32}
const _empty_lattice = _Lattice((NaN32, NaN32, NaN32, NaN32, NaN32, NaN32))
_is_complete(lattice::_Lattice) = !any(isnan, lattice)
_parse_lattice_value(s::AbstractString) = something(tryparse(Float32, s), NaN32)

# Lattice parameters from the CRYST1 field (PDB) or the _cell fields (mmCIF) of a file
function _read_lattice_parameters(io::IO)
    lattice = _empty_lattice
    for line in eachline(io)
        lattice = _read_lattice_line(lattice, line)
        _is_complete(lattice) && break
    end
    return lattice
end

# Updates the lattice parameters with the data of a CRYST1 (PDB) or _cell (mmCIF) line
function _read_lattice_line(lattice::_Lattice, line::AbstractString)
    a, b, c, α, β, γ = lattice
    fields = split(line)
    length(fields) >= 2 || return lattice
    # PDB file
    if fields[1] == "CRYST1" && length(fields) >= 7
        a = _parse_lattice_value(fields[2])
        b = _parse_lattice_value(fields[3])
        c = _parse_lattice_value(fields[4])
        α = _parse_lattice_value(fields[5])
        β = _parse_lattice_value(fields[6])
        γ = _parse_lattice_value(fields[7])
    end
    # mmCIF file
    if startswith(fields[1], "_cell.")
        val = @view(fields[1][7:end])
        val == "length_a" && (a = _parse_lattice_value(fields[2]))
        val == "length_b" && (b = _parse_lattice_value(fields[2]))
        val == "length_c" && (c = _parse_lattice_value(fields[2]))
        val == "angle_alpha" && (α = _parse_lattice_value(fields[2]))
        val == "angle_beta" && (β = _parse_lattice_value(fields[2]))
        val == "angle_gamma" && (γ = _parse_lattice_value(fields[2]))
    end
    return _Lattice((a, b, c, α, β, γ))
end

# Unit cell of a structure read from a file: `nothing` if the lattice parameters are
# incomplete, or if they define the 1Å unit cube, which is used in PDB files to indicate that 
# the structure is not crystallographic (e.g. NMR or cryo-EM structures).
function _structure_unitcell(lattice::_Lattice)
    _is_complete(lattice) || return nothing
    lattice == (1, 1, 1, 90, 90, 90) && return nothing
    return UnitCell(lattice_to_matrix(lattice...))
end


"""
    lattice_to_matrix(a, b, c, α, β, γ)

Converts unit cell lattice parameters and convert them to a 3x3
unit cell matrix (orthogonalization matrix), where the lattice vectors
are the *columns* of the matrix.

The resulting matrix has the lattice vectors as its columns, with
vector 'a' aligned along the x-axis and vector 'b' in the xy-plane.
This matrix can be used to transform fractional coordinates to
Cartesian coordinates.

# Arguments
- `a::Real`: Length of side a in Ångströms.
- `b::Real`: Length of side b in Ångströms.
- `c::Real`: Length of side c in Ångströms.
- `α::Real`: Angle alpha in degrees.
- `β::Real`: Angle beta in degrees.
- `γ::Real`: Angle gamma in degrees.

# Returns
- A 3x3 static matrix representing the unit cell vectors.

"""
function lattice_to_matrix(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real)
    a, b, c, α, β, γ = promote(float(a), float(b), float(c), float(α), float(β), float(γ))

    # Trigonometric functions of angles in degrees (exact for right angles)
    cos_α, cos_β, cos_γ = cosd(α), cosd(β), cosd(γ)
    sin_γ = sind(γ)

    # Calculate the unit cell volume V = abc * ω
    # ω is a term related to the Gram determinant
    ω = sqrt(1 - cos_α^2 - cos_β^2 - cos_γ^2 + 2 * cos_α * cos_β * cos_γ)

    # Construct the transformation matrix M
    T = typeof(a)
    M = SMatrix{3,3,T,9}(
        a, 
        zero(T),
        zero(T),
        b * cos_γ,
        b * sin_γ,
        zero(T),
        c * cos_β,
        c * (cos_α - cos_β * cos_γ) / sin_γ,
        c * ω / sin_γ
    )

    return M
end

"""
    matrix_to_lattice(M::AbstractMatrix)

Converts a 3x3 unit cell matrix where the lattice vectors are the *columns*, back into the six 
lattice parameters (sides a, b, c and angles α, β, γ).

# Arguments
- `M::AbstractMatrix`: A 3x3 matrix where columns are the lattice vectors.

# Returns
- `NamedTuple`: A named tuple containing the six parameters:
  `(a, b, c, α, β, γ)`. Lengths are in the matrix's original units
  (typically Ångströms), and angles are in degrees.

"""
function matrix_to_lattice(M::AbstractMatrix)
    T = eltype(M)
    # Ensure the input is a 3x3 matrix
    if size(M) != (3, 3)
        throw(ArgumentError("Input must be a 3x3 matrix."))
    end

    # Extract column vectors
    a_vec = @view(M[:, 1])
    b_vec = @view(M[:, 2])
    c_vec = @view(M[:, 3])

    # Calculate side lengths
    a = norm(a_vec)
    b = norm(b_vec)
    c = norm(c_vec)

    # Calculate angles in degrees
    # We use clamp() to avoid domain errors with acos due to floating point inaccuracies
    cos_gamma = clamp(dot(a_vec, b_vec) / (a * b), -one(T), one(T))
    cos_beta  = clamp(dot(a_vec, c_vec) / (a * c), -one(T), one(T))
    cos_alpha = clamp(dot(b_vec, c_vec) / (b * c), -one(T), one(T))

    α = rad2deg(acos(cos_alpha))
    β = rad2deg(acos(cos_beta))
    γ = rad2deg(acos(cos_gamma))

    return (a=a, b=b, c=c, α=α, β=β, γ=γ)
end

@testitem "read_unitcell" begin
    using PDBToolsBase

    a, b, c = 54.530, 63.780, 82.340
    α, β, γ = 90.0, 90.0, 90.0
    ortho_matrix = lattice_to_matrix(a, b, c, α, β, γ)
    @test ortho_matrix ≈ [54.530 0 0; 0 63.780 0; 0 0 82.340]
    @test eltype(ortho_matrix) == Float64
    @test all(v1 ≈ v2 for (v1,v2) in zip(values(matrix_to_lattice(ortho_matrix)),values((a = 54.53, b = 63.78, c = 82.34, α = 90.0, β = 90.0, γ = 90.0))))
    a, b, c = 54.530f0, 63.780f0, 82.340f0
    α, β, γ = 90.0f0, 90.0f0, 90.0f0
    ortho_matrix = lattice_to_matrix(a, b, c, α, β, γ)
    @test eltype(ortho_matrix) == Float32

    a, b, c = 39.850, 54.450, 59.840
    α, β, γ = 72.01, 89.44, 89.43
    triclinic_matrix = lattice_to_matrix(a, b, c, α, β, γ)
    @test triclinic_matrix ≈ [39.85 0.54168 0.584858; 0 54.4473 18.4767; 0 0 56.913] atol=1e-4
    @test eltype(triclinic_matrix) == Float64
    @test all(isapprox(v1,v2,atol=1e-5) for (v1,v2) in zip(values(matrix_to_lattice(triclinic_matrix)),values((a = 39.85, b = 54.449999999999996, c = 59.84, α = 72.01, β = 89.44, γ = 89.43))))

    @test read_unitcell(PDBToolsBase.TESTPBC) ≈ [107.845 0 0;  0.0  107.845 0; 0.0 0.0 107.845]
    @test read_unitcell(PDBToolsBase.TESTCIF) ≈ [74.449 -3.25427f-6 -3.25427f-6; 0.0 74.449 -3.25427f-6; 0.0 0.0 74.449]
    
    # Throw error if unitcel ldata is missing
    @test_throws ArgumentError read_unitcell(PDBToolsBase.TESTPDB)
    @test_throws ArgumentError read_unitcell(PDBToolsBase.BROKENCIF)
    @test_throws ArgumentError matrix_to_lattice([1 0; 0 1])

    # Test allocations of conversions
    using Chairmarks
    a, b, c = 39.850, 54.450, 59.840
    α, β, γ = 72.01, 89.44, 89.43
    t = @b lattice_to_matrix($a, $b, $c, $α, $β, $γ)
    @test t.allocs == 0
    l = lattice_to_matrix(a, b, c, α, β, γ)
    t = @b matrix_to_lattice($l)
    @test t.allocs == 0

end