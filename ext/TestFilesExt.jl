module TestFilesExt

using LazyArtifacts: LazyArtifacts, @artifact_str
import PDBToolsBase: test_file

test_file(name::AbstractString) = joinpath(artifact"pdbtools_test_files", name)

end
