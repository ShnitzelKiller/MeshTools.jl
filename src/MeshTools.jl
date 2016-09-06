module MeshTools
export convexhull, createMesh, findFacet, removeDoubles, Node, saveObj, volume, surfacearea, MeshGrid, countParts, separate, linearFilter, makeGaussianFilter, joinmeshes, Mesh, cylinder, invertnormals, saveOff

using Compat
import Compat.String
typealias Mesh Tuple{Matrix{Float64}, Vector{Int}}
include("meshfit.jl")
include("meshvolume.jl")
include("hull.jl")
include("objconverter.jl")
include("meshops.jl")
include("primitives.jl")
end
