module MeshTools
export convexhull, createMesh, findFacet, removeDoubles, Node, saveObj, volume, surfacearea, MeshGrid, countParts, separate, linearFilter, makeGaussianFilter
include("meshfit.jl")
include("meshvolume.jl")
include("hull.jl")
include("objconverter.jl")
end
