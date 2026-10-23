using TestItemRunner

@testmodule AllocTest begin
    # This module defines the Allocs struct and the comparison operators
    # to conditionally compare the number of allocations based on the
    # BUILD_IS_PRODUCTION_BUILD environment variable.
    export Allocs
    @kwdef struct Allocs
        prodbuild::Bool = haskey(ENV, "BUILD_IS_PRODUCTION_BUILD") && ENV["BUILD_IS_PRODUCTION_BUILD"] == "true"
        allocs::Int
    end
    Allocs(allocs::Int) = Allocs(; allocs)
    import Base: ==, >, <
    ==(a::Real, b::Allocs) = b.prodbuild ? Int(a) == b.allocs : true
    <(a::Real, b::Allocs) = b.prodbuild ? Int(a) < b.allocs : true
    ==(a::Allocs, b::Real) = a.prodbuild ? a.allocs == Int(b) : true
    <(a::Allocs, b::Real) = a.prodbuild ? a.allocs < Int(b) : true
end

# PDBToolsBase exports nothing: this snippet brings its public names into the scope
# of the test items
@testsnippet PublicNames begin
    using LazyArtifacts # loads the extension that defines test_file
    for n in names(PDBToolsBase)
        n == :PDBToolsBase || @eval using PDBToolsBase: $n
    end
end

@run_package_tests

@testitem "Aqua.test_all" setup = [PublicNames] begin
    import Aqua
    Aqua.test_all(PDBToolsBase)
end

# The public API is defined in docs/src/api.md: every public name must have a docstring
# (the documentation build checks that all of them are included in the API page)
@testitem "public names are documented" setup = [PublicNames] begin
    public_names = filter(n -> Base.ispublic(PDBToolsBase, n), names(PDBToolsBase; all=true))
    @test length(public_names) > 50
    @test filter(n -> !Base.Docs.hasdoc(PDBToolsBase, n), public_names) == Symbol[]
end
