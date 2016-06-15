# MeshTools
This package provides a variety of mesh and volume-related tools, including
- Marching cubes triangulation
- Convex hull generation
- Mesh volume computation
The meshes used by functions in this package are in the index array format. A mesh is therefore a tuple, where the first element is a 3xn array whose columns represent vertices, and the second element is a list of integer indices, where each group of three indices looks up a triangle in the vertex array.