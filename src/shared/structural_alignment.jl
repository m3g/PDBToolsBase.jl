"""
    center_of_mass(x::AbstractVector{<:AbstractVector}[, mass::AbstractVector=nothing])

Calculate the center of mass of a set of points.

# Arguments

- `x::AbstractVector{<:AbstractVector}`: A vector of coordinates.
- `mass::AbstractVector`: A vector of masses. If not provided, all masses are assumed to be equal.

# Example

```jldoctest
julia> import PDBTools: center_of_mass

julia> x = [ [1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0] ];

julia> center_of_mass(x)
3-element Vector{Float64}:
 4.0
 5.0
 6.0

julia> center_of_mass(x, [1.0, 2.0, 3.0]) # providing masses
3-element Vector{Float64}:
 5.0
 6.0
 7.0

```

"""
center_of_mass(x::AbstractVector{<:AbstractVector}) = center_of_mass(x, nothing)
center_of_mass(x::AbstractVector{<:AbstractVector}, ::Nothing) = sum(x) / length(x)
center_of_mass(x::AbstractVector{<:AbstractVector}, mass::AbstractVector) =
    sum(x[i] * mass[i] for i in eachindex(x, mass)) / sum(mass)

"""
    rmsd(x::AbstractVector,y::AbstractVector)

Calculate the root mean square deviation between two vectors of coordinates.

# Arguments

- `x::AbstractVector`: A vector of coordinates.
- `y::AbstractVector`: A vector of coordinates.

# Returns

- `rmsd::Real`: The root mean square deviation between the two vectors.

```jldoctest
julia> import PDBTools: rmsd

julia> x = [ [1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0] ];

julia> y = [ [2.0, 3.0, 4.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0] ];

julia> rmsd(x, y)
1.0
```

"""
function rmsd(x::AbstractVector{<:AbstractVector}, y::AbstractVector{<:AbstractVector})
    rmsd = zero(eltype(first(x)))
    for i in eachindex(x, y)
        rmsd += sum(abs2, x[i] .- y[i])
    end
    return sqrt(rmsd / length(x))
end

"""
    align(x, y; mass=nothing)
    align!(x, y; mass=nothing)

Aligns two structures (sets of points in 3D space). Solves the "Procrustes" problem, which is to find the best
translation, and rotation, that aligns the two structures, minimizing the RMSD between them.

Structures are expected to be of the same size, and the correspondence is assumed from the vector indices. 

`align` returns a new vector containing the coordinates of x aligned to y. 
`align!` modifies the input vector `x` in place.

"""
function align end
@doc (@doc align) align!

function align(
    x::AbstractVector{<:AbstractVector},
    y::AbstractVector{<:AbstractVector};
    mass=nothing
)
    xnew = copy(x)
    return align!(xnew, y; mass)
end

"""
    apply_alignment_transformation!(x, cmx, cmy, rotation_matrix)

Translates and rotates `x` such that the coordinates `cmx` are moved to `cmy`, and
such apply the rotation defined in `rotation_matrix`. 

"""
function apply_alignment_transformation!(x, cmx, cmy, rotation_matrix)
    for i in eachindex(x)
        x[i] = rotation_matrix * (x[i] - cmx) + cmy
    end
    return x
end

"""
    alignment_movements(
        x::AbstractVector{<:AbstractVector},
        y::AbstractVector{<:AbstractVector};
        mass=nothing,
        xm=zeros(3, length(x)),
        xp=zeros(3, length(x)),
    )
    
Returns `cmx`, `cmy`, and `rotation_matrix`, necessary to translate and rotate y to minimize the RMSD relative to x.

Can be used to obtain the movement that aligns two sets of points, of the same size, and apply these movements to 
another set of points.

"""
function alignment_movements(
    x::AbstractVector{<:AbstractVector},
    y::AbstractVector{<:AbstractVector};
    mass=nothing,
    xm=zeros(3, length(x)),
    xp=zeros(3, length(x)),
)
    length(x) == length(y) || throw(DimensionMismatch("x and y must have the same length"))
    (length(x[1]) != 3 || length(x[2]) != 3) && throw(DimensionMismatch("x and y must be 3D vectors"))

    cmx = center_of_mass(x, mass)
    cmy = center_of_mass(y, mass)
    for i in eachindex(x, y)
        x[i] -= cmx
        y[i] -= cmy
    end

    for i in eachindex(x, y)
        xm[1:3, i] .= y[i] .- x[i]
        xp[1:3, i] .= y[i] .+ x[i]
    end

    q = zeros(MMatrix{4,4,eltype(xm),16})
    for i in eachindex(x)
        q[1, 1] = q[1, 1] + sum(abs2, @view(xm[1:3, i]))
        q[1, 2] = q[1, 2] + xp[2, i] * xm[3, i] - xm[2, i] * xp[3, i]
        q[1, 3] = q[1, 3] + xm[1, i] * xp[3, i] - xp[1, i] * xm[3, i]
        q[1, 4] = q[1, 4] + xp[1, i] * xm[2, i] - xm[1, i] * xp[2, i]
        q[2, 2] = q[2, 2] + xp[2, i]^2 + xp[3, i]^2 + xm[1, i]^2
        q[2, 3] = q[2, 3] + xm[1, i] * xm[2, i] - xp[1, i] * xp[2, i]
        q[2, 4] = q[2, 4] + xm[1, i] * xm[3, i] - xp[1, i] * xp[3, i]
        q[3, 3] = q[3, 3] + xp[1, i]^2 + xp[3, i]^2 + xm[2, i]^2
        q[3, 4] = q[3, 4] + xm[2, i] * xm[3, i] - xp[2, i] * xp[3, i]
        q[4, 4] = q[4, 4] + xp[1, i]^2 + xp[2, i]^2 + xm[3, i]^2
    end
    q[2, 1] = q[1, 2]
    q[3, 1] = q[1, 3]
    q[3, 2] = q[2, 3]
    q[4, 1] = q[1, 4]
    q[4, 2] = q[2, 4]
    q[4, 3] = q[3, 4]
    q = SMatrix(q)

    # Computing the eigenvectors 'v' of the q matrix
    v = eigvecs(q)

    # Compute rotation matrix
    u = zeros(MMatrix{3,3,Float64,9})
    u[1, 1] = v[1, 1]^2 + v[2, 1]^2 - v[3, 1]^2 - v[4, 1]^2
    u[1, 2] = 2.0 * (v[2, 1] * v[3, 1] + v[1, 1] * v[4, 1])
    u[1, 3] = 2.0 * (v[2, 1] * v[4, 1] - v[1, 1] * v[3, 1])
    u[2, 1] = 2.0 * (v[2, 1] * v[3, 1] - v[1, 1] * v[4, 1])
    u[2, 2] = v[1, 1]^2 + v[3, 1]^2 - v[2, 1]^2 - v[4, 1]^2
    u[2, 3] = 2.0 * (v[3, 1] * v[4, 1] + v[1, 1] * v[2, 1])
    u[3, 1] = 2.0 * (v[2, 1] * v[4, 1] + v[1, 1] * v[3, 1])
    u[3, 2] = 2.0 * (v[3, 1] * v[4, 1] - v[1, 1] * v[2, 1])
    u[3, 3] = v[1, 1]^2 + v[4, 1]^2 - v[2, 1]^2 - v[3, 1]^2
    u = SMatrix(u)

    # Restore arrays
    for i in eachindex(x, y)
        x[i] += cmx
        y[i] += cmy
    end

    return cmx, cmy, u
end

function align!(
    x::AbstractVector{<:AbstractVector},
    y::AbstractVector{<:AbstractVector};
    mass=nothing,
    # Auxiliary arrays that might be preallocated
    xm=zeros(3, length(x)),
    xp=zeros(3, length(x))
)
    cmx, cmy, u = alignment_movements(x, y; mass, xm, xp)
    apply_alignment_transformation!(x, cmx, cmy, u)
    return x
end

@testitem "structural_alignment" setup = [PublicNames] begin
    import PDBToolsBase: center_of_mass, align, rmsd
    using StaticArrays: SVector
    using Rotations: RotMatrix3

    x = [rand(SVector{3,Float64}) for _ in 1:10]
    @test center_of_mass(x) ≈ sum(x) / length(x)

    y = x .+ Ref(SVector{3}(1, 1, 1))
    @test rmsd(x, y) ≈ sqrt(length(x) * 3 / length(x))

    # apply a random rotation and translation to x
    y = x .+ Ref(SVector{3}(45.0, -15.0, 31.5))
    y .= Ref(rand(RotMatrix3)) .* y
    @test rmsd(x, y) > 0.0
    z = align(x, y)
    @test rmsd(z, y) ≈ 0.0 atol = 1e-5

end


