# MeshTools
```jl
Pkg.clone("https://github.com/ShnitzelKiller/MeshTools.jl.git")
```

This package provides a variety of mesh and volume-related tools, including
- Marching cubes triangulation
- Convex hull generation
- Mesh volume computation
- Scalar field generation from weighted points

The meshes used by functions in this package are in the index array format. A mesh is therefore a tuple, where the first element is a 3xn array whose columns represent vertices, and the second element is a list of integer indices, where each group of three indices looks up a triangle in the vertex array.

## Examples
Here are examples of how to work with meshes. Lists of points are given in Matlab-like format, as 3xn matrices. If you just have a collection of points, the below conversion is necessary for use with this package:
```jl
julia> points = hcat(points...)
```
Now we can create a convex hull using `points`:
```jl
julia> verts, inds = convexhull(points)
```

Other methods, such as the `MeshGrid` constructor, require weighted points, which are represented as a 4xn matrix. To create such a matrix from a list of points and a list of weights, run
```jl
julia> points = hcat(points...)
julia> wpoints = vcat(points, energies')
```
now `wpoints` can be used in the constructor:
```jl
julia> grid = MeshGrid(25, 25, 25, wpoints)
```

## Summary
### MeshGrid
`MeshGrid(resX, resY, resZ, hits, dilate=false, stdev=20.0, buffer=[0.0, 0.0, 0.0], threshold=Inf)` takes a set of weighted points and bins them into a 3D array. It optionally applies a Gaussian blur to the data, making it ideal for use with `createMesh`. If this is the goal, supplying the threshold used for marching cubes triangulation to the constructor allows it to resize the volume boundaries to ensure that the mesh is fully contained. The `buffer` parameter is just an initial guess for how much extra space you will need, as a percentage of the sidelength in each dimension. If `threshold` is set, `buffer` will always be recalculated if more space is needed.
The `MeshGrid` object can be indexed by real position:
```jl
julia> grid[0.3, 5.2, 9.9]
8.8
```
Values are linearly interpolated so as to give continuous results.
### Convex Hull
`convexhull(points)` creates the convex hull of the provided points matrix as a triangle mesh.
### Marching Cubes
`createMesh{T}(data :: Array{T, 3}, scaleX, scaleY, scaleZ, originX, originY, originZ, threshold, startIndex=1, uniformScale=1)` creates the marching cubes triangulation of the 3D array `data` with the given dimensions and origin point, with an optional shifted mesh index position and uniform scale. The former might be useful when creating mesh parts with the intention of merging them together. A simplified constructor exists for MeshGrid objects which already contain position information:
`createMesh(grid :: MeshGrid, threshold, startIndex=1, uniformScale=1)`
### Analysis Tools
A number of functions are provided which act on meshes:
```jl
volume(verts, inds)
surfacearea(verts, inds)
verts, inds = joinmeshes(meshes...) #where each mesh is a tuple of verts and indices
verts, inds = invertnormals(verts, inds)
```
