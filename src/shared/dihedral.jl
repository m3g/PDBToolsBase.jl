"""
    dihedral(v1, v2, v3, v4; degrees=true)
    dihedral(v::AbstractVector; degrees=true)

Computes the dihedral angle between the planes formed by the vectors `v1-v2` and `v2-v3`, and `v2-v3` and `v3-v4`.
The input vectors must have 3 elements. The function returns the dihedral angle in radians or degrees.

If the input is a vector with 4 vectors, the function computes the dihedral angle between the planes 
formed by the vectors `v[1]-v[2]` and `v[2]-v[3]`, and `v[2]-v[3]` and `v[3]-v[4]`.

The optional argument `degrees` specifies whether the output is in degrees (default) or radians.

"""
function dihedral(
    v1::AbstractVector, 
    v2::AbstractVector, 
    v3::AbstractVector, 
    v4::AbstractVector; 
    degrees=true
)
    if any(length(v) != 3 for v in (v1,v2,v3,v4))
        throw(ArgumentError("All input vectors must have 3 elements"))
    end 

    # Vectors connecting consecutive atoms
    v21 = v2 - v1 
    v32 = v3 - v2 
    v43 = v4 - v3

    # Normal vectors to the planes formed by the vectors v1 and v2, and v2 and v3
    n1 = cross(v21, v32)
    n2 = cross(v32, v43)

    n1 = n1 / norm(n1)
    n2 = n2 / norm(n2)

    # Normalized vector perpendicular to the plane formed by the vectors n1 and n2
    u1 = cross(n1, n2) / norm(n1)

    # Normalize v2
    unitary_v32 = v32 / norm(v32)

    # Computes the projection of the vector u1 in the direction of unitary_v2,
    # and of n1 in the direction of n2
    m1 = dot(u1, unitary_v32)  
    m2 = dot(n1, n2)         

    # Computes the dihedral angle in radians or degrees
    return degrees ? atand(m1, m2) : atan(m1, m2)
end

function dihedral(v::AbstractVector; degrees=true)
    length(v) == 4 || throw(ArgumentError("The input vector or vectors must have 4 elements"))
    return dihedral(v[1], v[2], v[3], v[4]; degrees)
end

@testitem "dihedral" begin
    using PDBToolsBase: dihedral
    using StaticArrays

    d = dihedral( 
        Float32[-9.229, -14.861, -5.481], 
        Float32[-10.048, -15.427, -5.569], 
        Float32[-9.488, -13.913, -5.295], 
        Float32[-8.652, -15.208, -4.741],
    )
    @test d ≈ -34.57 atol=1e-2

    d = dihedral(
        Float32[-8.483, -14.912, -6.726], 
        Float32[-5.113, -13.737, -5.466], 
        Float32[-3.903, -11.262, -8.062], 
        Float32[-1.162, -9.64, -6.015],
    )
    @test d ≈ 164.43 atol=1e-2

    d = dihedral(
        Float32[-9.229, -14.861, -5.481], 
        Float32[-8.483, -14.912, -6.726], 
        Float32[-7.227, -14.047, -6.599], 
        Float32[-7.083, -13.048, -7.303],
    )
    @test d ≈ -115.83 atol=1e-2

    d2 = dihedral(
        Float32[-9.229, -14.861, -5.481], 
        Float32[-8.483, -14.912, -6.726], 
        Float32[-7.227, -14.047, -6.599], 
        Float32[-7.083, -13.048, -7.303];
        degrees=false
    )
    @test d2 ≈ deg2rad(d) 

    @test dihedral(
        SVector(-9.229, -14.861, -5.481), 
        SVector(-8.483, -14.912, -6.726), 
        SVector(-7.227, -14.047, -6.599), 
        SVector(-7.083, -13.048, -7.303);
    ) ≈ d

end

"""
    dihedrals(v::AbstractVector{<:AbstractVector}; degrees=true)

Computes the dihedral angles for many sets of 4 vectors. The input is a vector of vectors, where each
element is a vector with 4 vectors. The function returns a vector with the dihedral angles in radians or degrees.

## Example

```jldoctest; filter = r"([0-9]+\\.[0-9]{2})[0-9]+" => s"\\1***"
julia> using PDBTools: dihedrals

julia> v1 = [[-8.483, -14.912, -6.726], [-5.113, -13.737, -5.466], [-3.903, -11.262, -8.062], [-1.162, -9.64, -6.015]];

julia> v2 = [[-9.229, -14.861, -5.481], [-8.483, -14.912, -6.726], [-7.227, -14.047, -6.599], [-7.083, -13.048, -7.303]];

julia> dihedrals([v1,v2])
2-element Vector{Float64}:
  164.43481280739516
 -115.82544005374316
```

"""
function dihedrals(v::AbstractVector{<:AbstractVector}; degrees=true)
    ds = zeros(eltype(eltype(eltype(v))), length(v))
    for (i, v4) in enumerate(v)
        ds[i] += dihedral(v4; degrees)
    end
    return ds
end

@testitem "dihedrals" begin
    using PDBToolsBase: dihedral, dihedrals
    using StaticArrays
    x = SVector{3, Float32}[[-8.483, -14.912, -6.726], [-5.113, -13.737, -5.466], [-3.903, -11.262, -8.062], [-1.162, -9.64, -6.015]]
    d = dihedral(x)
    x10 = [ x for _ in 1:10 ]
    ds = dihedrals(x10)
    @test all(ds .≈ d)
end
