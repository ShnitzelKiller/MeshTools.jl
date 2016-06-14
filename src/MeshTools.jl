module MeshTools
export convexhull, createMesh, findFacet, removeDoubles, Node, saveObj, volume, MeshGrid
include("meshfit.jl")
include("meshvolume.jl")
include("meshgrid.jl")
include("hull.jl")
include("objconverter.jl")
end
