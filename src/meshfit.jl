include("octree.jl")
include("meshgrid.jl")

using DataStructures

typealias Voxel{T} Array{T, 3}

const v0 = 0x01 :: UInt8
const v1 = 0x02 :: UInt8
const v2 = 0x04 :: UInt8
const v3 = 0x08 :: UInt8
const v4 = 0x10 :: UInt8
const v5 = 0x20 :: UInt8
const v6 = 0x40 :: UInt8
const v7 = 0x80 :: UInt8

#The ith columns is the two vertex indices adjacent to the ith edge.
const edgeCouplings =
[
0 1 2 0 4 5 6 4 0 1 2 3;
1 2 3 3 5 6 7 7 4 5 6 7
]

#If the 8-bit bitmask representing the covered vertices is used as an index,
#each 12-bit value represents the combination of edges intersected by the isosurface.
const edgeTable = [
0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   ]

#If the 8-bit bitmask representing the covered vertices is used as an index,
#the resulting column is the sequence of edges visited by the triangles constructing the isosurface.
const triTable =
[
0 1 1 2 2 1 10 3 4 1 2 2 4 1 4 10 5 5 1 5 2 4 10 3 9 12 10 5 4 2 5 5 10 10 1 9 2 4 6 3 10 1 1 3 11 5 6 6 10 10 1 2 10 11 9 3 8 10 3 12 10 6 12 12 11 1 10 2 2 2 10 6 3 12 1 6 7 1 4 7 6 5 2 11 7 2 9 8 4 6 1 10 9 6 1 7 11 5 11 9 2 4 1 9 11 1 4 7 10 9 4 7 8 1 11 11 2 3 8 8 3 3 2 12 9 1 8 8 8 4 1 9 11 2 3 7 8 8 3 2 11 11 1 8 7 4 9 10 7 2 5 11 9 1 2 2 9 11 5 11 5 1 6 12 10 7 8 4 8 10 4 7 10 2 5 8 7 4 1 7 2 1 12 7 6 10 2 2 2 11 1 11 12 12 6 11 12 1 10 8 3 9 10 10 2 1 10 10 6 6 1 11 3 1 1 10 3 6 4 6 9 1 9 10 5 1 2 4 5 10 12 12 3 10 4 2 5 5 5 5 10 4 1 4 2 4 1 4 3 10 3 2 2 1 1 0;
0 9 2 9 3 9 3 9 12 12 10 12 11 11 10 9 8 4 2 2 3 5 3 11 5 5 1 8 11 12 8 8 6 6 6 6 3 1 3 11 6 12 6 2 4 10 5 5 8 4 8 6 8 2 1 11 10 6 4 3 6 8 11 11 7 9 1 9 7 7 7 10 4 1 2 11 4 9 12 6 11 4 10 7 2 3 5 4 12 11 2 3 5 2 6 6 5 11 1 4 5 1 3 4 5 9 12 5 7 12 12 5 11 8 7 7 3 7 9 4 4 1 9 3 10 10 9 12 7 1 2 2 2 3 10 12 3 1 8 7 8 8 4 7 9 7 7 5 9 3 12 10 3 5 10 10 2 2 7 10 10 9 1 8 6 12 7 5 3 6 7 3 6 7 1 7 10 7 12 12 3 12 9 12 9 6 6 6 4 2 4 6 6 6 12 8 2 9 8 6 6 3 1 9 4 9 1 9 9 1 2 12 6 5 3 5 6 11 11 11 5 5 5 5 12 9 11 2 12 8 8 8 10 11 8 11 10 10 1 9 11 1 2 2 3 1 3 3 4 11 4 11 4 10 4 0;
0 4 10 4 11 4 11 4 3 3 1 3 2 2 1 11 9 1 10 10 11 8 11 10 8 8 2 12 2 11 9 12 5 5 5 5 11 9 11 6 5 3 5 6 12 6 1 9 9 1 9 4 9 3 3 6 6 8 12 2 9 1 1 6 6 4 2 4 6 6 6 9 12 9 10 7 12 12 7 10 7 1 1 6 3 6 8 10 3 7 10 2 8 12 10 10 10 7 2 2 10 9 5 3 10 3 3 2 5 2 7 9 7 4 8 8 7 10 1 3 12 8 1 2 7 2 1 7 12 9 10 10 3 11 1 8 4 9 7 3 7 7 8 11 5 12 12 7 5 11 9 4 4 3 1 5 4 1 4 5 6 4 2 7 5 8 12 9 4 5 3 9 5 11 11 11 6 12 9 4 11 4 6 4 10 7 9 7 7 1 9 7 11 11 8 6 3 4 6 3 11 1 2 3 6 8 4 8 5 5 10 5 2 12 6 6 11 3 3 3 6 6 6 6 8 4 12 5 8 5 5 5 11 8 11 3 2 2 4 8 9 10 11 11 12 10 12 12 9 3 9 3 9 2 9 0;
0 0 0 10 0 2 1 3 0 9 3 2 12 1 4 11 0 8 9 5 9 4 10 3 4 12 9 10 4 2 10 5 0 1 2 9 10 2 6 4 3 1 1 3 11 1 6 6 6 10 1 4 10 10 9 3 8 10 1 12 9 6 12 8 0 6 6 2 3 2 10 6 11 12 3 2 7 1 1 7 5 5 6 2 7 6 10 8 8 5 5 10 4 6 1 7 7 5 11 9 2 2 5 9 11 3 1 7 10 9 4 12 8 1 2 11 2 3 8 7 11 3 2 12 9 12 8 0 0 12 12 9 7 4 3 3 7 8 3 2 11 2 1 8 12 4 9 10 7 4 5 11 9 5 3 2 9 11 5 7 8 5 6 9 11 2 6 4 8 1 4 7 11 2 5 8 7 1 1 7 10 1 12 7 6 10 2 3 2 11 6 0 8 12 6 11 12 2 10 8 3 9 6 10 4 1 10 6 6 6 9 11 3 1 1 3 4 6 4 6 9 2 9 0 5 5 2 4 10 10 12 12 3 10 4 9 5 5 8 0 11 4 1 12 2 4 9 0 3 1 3 0 10 0 0 0;
0 0 0 9 0 3 3 11 0 12 4 10 11 9 12 9 0 4 5 8 5 1 1 10 12 3 5 5 12 5 1 12 0 9 6 4 6 3 5 3 4 9 2 6 2 9 1 9 8 6 2 6 6 6 3 6 9 8 2 2 6 1 1 12 0 11 11 10 7 3 1 9 7 3 4 10 6 12 4 10 8 8 11 10 6 3 1 10 9 8 8 12 12 12 7 10 5 10 7 2 3 3 3 3 7 9 2 2 4 2 7 7 9 11 11 8 7 10 1 8 7 8 8 2 7 7 1 0 0 8 8 4 12 1 11 11 3 7 4 9 2 8 8 11 9 1 5 7 12 1 7 4 5 7 4 5 7 1 4 11 7 10 5 4 2 3 5 6 7 9 8 9 2 8 11 11 12 7 6 4 6 7 6 4 3 7 9 2 7 1 7 0 6 8 11 12 8 3 3 3 4 6 11 3 8 8 4 10 11 12 5 5 9 12 6 12 6 3 6 3 6 1 6 0 10 10 12 5 12 12 5 5 8 8 11 8 2 2 5 0 12 10 11 4 12 10 1 0 9 10 9 0 2 0 0 0;
0 0 0 2 0 11 10 9 0 1 12 12 4 11 10 12 0 5 8 2 8 5 3 8 3 5 8 12 11 12 12 10 0 4 1 6 5 11 3 6 12 12 6 9 4 2 12 11 10 4 8 8 8 1 6 4 10 3 9 8 8 10 4 6 0 7 7 9 2 7 7 3 6 1 12 3 4 6 7 12 9 4 7 8 2 7 6 5 5 3 9 3 6 7 6 12 11 11 1 7 5 10 7 5 5 12 7 11 7 1 1 9 11 8 8 2 9 2 7 3 9 12 9 8 8 8 7 0 0 7 7 2 8 9 10 4 8 1 8 7 8 11 11 9 7 7 7 4 9 12 12 3 3 3 5 3 2 7 9 5 12 6 1 5 3 11 11 5 3 7 7 8 7 7 6 9 10 4 12 6 12 12 7 6 9 1 1 7 11 7 11 0 12 6 12 8 2 8 8 12 6 3 4 2 6 2 6 8 9 1 11 6 6 4 10 4 3 5 11 5 4 6 4 0 12 8 5 9 5 8 3 3 10 5 3 5 8 8 4 0 9 12 9 11 10 12 12 0 11 3 11 0 9 0 0 0;
0 0 0 0 0 0 0 11 0 0 0 10 0 9 12 0 0 0 0 8 0 2 9 3 0 3 3 10 8 2 10 10 0 0 0 4 0 5 5 4 0 5 3 3 10 9 6 11 0 6 2 0 11 6 9 4 4 10 2 8 11 8 11 0 0 0 0 6 0 4 1 6 0 11 6 10 6 1 1 12 0 7 9 2 5 4 1 4 11 5 3 10 4 2 1 5 0 1 7 9 3 3 0 5 12 5 1 5 10 12 1 0 9 1 2 2 2 7 7 0 11 1 2 11 10 0 4 0 0 0 0 12 0 7 7 11 0 7 1 2 2 2 1 9 0 1 10 10 3 1 1 10 5 0 3 3 9 7 7 0 0 12 8 4 8 1 5 4 6 1 2 3 2 2 1 6 12 1 1 6 10 1 9 3 6 1 6 0 4 10 0 0 0 9 2 10 8 2 10 6 4 9 6 9 0 2 6 0 11 6 9 12 3 5 3 0 4 5 4 2 4 0 10 0 10 10 2 2 10 10 3 9 3 11 8 0 8 1 0 0 0 12 9 0 10 2 0 0 11 0 1 0 0 0 0 0;
0 0 0 0 0 0 0 10 0 0 0 9 0 12 11 0 0 0 0 4 0 3 5 8 0 1 4 12 9 1 12 12 0 0 0 2 0 10 1 6 0 10 4 9 6 11 12 9 0 8 6 0 2 4 6 6 12 3 8 2 2 12 6 0 0 0 0 11 0 1 3 3 0 7 11 12 2 6 7 10 0 6 5 8 8 1 7 3 7 3 4 5 6 1 4 8 0 9 5 7 7 5 0 3 3 10 7 9 2 7 7 0 10 10 8 8 9 8 1 0 9 10 11 7 2 0 12 0 0 0 0 8 0 12 12 9 0 3 2 10 4 9 11 11 0 5 1 4 11 7 3 5 7 0 5 5 5 1 11 0 0 8 7 6 7 9 3 3 5 7 6 2 8 1 4 5 9 6 2 4 12 10 1 11 7 7 7 0 9 6 0 0 0 4 10 9 6 8 1 10 8 8 4 8 0 8 4 0 12 11 11 4 12 6 12 0 5 3 9 10 6 0 1 0 11 12 5 11 3 2 5 4 4 3 5 0 2 9 0 0 0 10 11 0 12 3 0 0 9 0 2 0 0 0 0 0;
0 0 0 0 0 0 0 9 0 0 0 12 0 11 10 0 0 0 0 2 0 11 8 4 0 5 12 3 5 5 11 11 0 0 0 6 0 6 3 5 0 6 12 12 5 2 11 12 0 4 8 0 3 1 8 8 3 1 9 6 4 1 1 0 0 0 0 7 0 9 7 7 0 6 7 3 4 2 6 9 0 11 8 4 9 5 6 10 6 1 12 12 2 12 7 10 0 4 1 5 5 10 0 7 4 11 5 2 4 2 5 0 11 11 9 4 10 10 3 0 10 8 8 2 7 0 1 0 0 0 0 7 0 8 8 4 0 1 10 9 8 8 10 10 0 7 2 2 2 12 10 4 3 0 7 7 7 5 4 0 0 7 12 5 12 4 11 6 10 3 1 9 7 8 11 11 10 7 6 2 9 7 6 4 3 3 9 0 7 1 0 0 0 1 1 2 2 6 3 3 6 6 8 3 0 6 8 0 9 12 12 5 9 12 6 0 6 1 6 3 2 0 6 0 12 8 1 5 12 12 1 5 8 8 11 0 4 2 0 0 0 11 12 0 9 10 0 0 10 0 9 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 0 0 0 10 0 8 12 0 0 0 0 0 0 0 0 4 0 0 0 5 0 9 12 0 0 0 0 0 0 6 11 0 0 3 2 0 11 2 9 0 0 0 0 0 0 0 0 4 0 0 0 10 0 6 1 0 0 0 0 8 0 4 1 6 0 3 6 8 6 8 12 8 0 0 0 7 0 3 0 0 0 5 7 3 12 10 0 0 0 7 2 0 9 1 0 0 9 7 7 7 12 0 12 0 0 0 0 0 0 0 0 11 0 0 0 9 0 2 7 0 0 0 0 12 0 1 3 12 0 0 5 0 7 0 1 0 0 0 0 4 0 5 5 11 0 7 6 5 2 9 7 5 0 1 6 0 12 6 11 11 4 0 4 0 6 6 0 0 0 0 0 9 0 8 3 4 0 11 4 11 0 0 0 0 0 12 11 10 5 3 5 0 4 0 5 10 0 0 1 0 0 10 8 8 10 3 0 4 8 9 2 0 0 9 0 0 0 0 0 0 0 3 0 0 0 0 2 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10 0 0 0 3 0 12 1 0 0 0 0 0 0 0 0 5 0 0 0 9 0 12 1 0 0 0 0 0 0 8 6 0 0 8 6 0 4 1 1 0 0 0 0 0 0 0 0 3 0 0 0 9 0 12 6 0 0 0 0 10 0 5 3 10 0 8 11 12 12 12 7 12 0 0 0 2 0 7 0 0 0 11 2 2 7 2 0 0 0 8 9 0 7 10 0 0 7 8 8 8 7 0 7 0 0 0 0 0 0 0 0 10 0 0 0 8 0 1 11 0 0 0 0 4 0 5 11 4 0 0 4 0 11 0 4 0 0 0 0 2 0 10 1 6 0 9 5 9 4 8 11 9 0 10 7 0 6 7 6 6 9 0 9 0 7 7 0 0 0 0 0 4 0 3 12 3 0 3 11 3 0 0 0 0 0 4 5 5 6 12 6 0 9 0 6 5 0 0 4 0 0 11 5 5 2 12 0 3 5 8 11 0 0 8 0 0 0 0 0 0 0 12 0 0 0 0 11 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 0 0 0 2 0 5 4 0 0 0 0 0 0 0 0 9 0 0 0 6 0 11 4 0 0 0 0 0 0 4 3 0 0 12 8 0 12 11 8 0 0 0 0 0 0 0 0 9 0 0 0 12 0 7 10 0 0 0 0 5 0 8 7 7 0 12 7 5 7 5 4 10 0 0 0 11 0 5 0 0 0 7 11 12 4 5 0 0 0 11 1 0 8 4 0 0 8 11 11 2 4 0 1 0 0 0 0 0 0 0 0 9 0 0 0 7 0 9 8 0 0 0 0 7 0 7 10 7 0 0 9 0 2 0 10 0 0 0 0 6 0 6 3 3 0 8 1 6 8 1 8 11 0 6 12 0 7 10 3 4 3 0 3 0 10 1 0 0 0 0 0 2 0 12 8 9 0 6 3 6 0 0 0 0 0 1 6 2 9 2 9 0 5 0 9 3 0 0 6 0 0 12 12 12 3 2 0 5 10 1 1 0 0 2 0 0 0 0 0 0 0 10 0 0 0 0 9 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 12 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 6 0 1 9 0 0 0 0 0 0 0 0 0 0 0 0 9 0 7 0 0 0 0 0 0 0 8 0 0 0 10 3 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 12 0 0 0 2 0 10 4 0 0 0 0 0 0 2 1 0 0 0 7 0 9 0 0 0 0 0 0 0 0 0 0 10 0 0 0 8 0 0 0 0 0 0 0 4 0 6 12 0 0 0 1 0 0 0 0 0 0 0 0 11 0 1 0 0 0 3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 0 0 0 11 0 5 5 0 0 0 0 0 0 0 0 0 0 0 0 12 0 5 0 0 0 0 0 0 0 4 0 0 0 11 4 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 0 0 0 0 0 0 10 0 0 0 0 0 0 0 0 8 0 0 0 6 0 6 8 0 0 0 0 0 0 3 3 0 0 0 3 0 10 0 0 0 0 0 0 0 0 0 0 9 0 0 0 6 0 0 0 0 0 0 0 2 0 2 9 0 0 0 2 0 0 0 0 0 0 0 0 12 0 9 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10 0 0 0 7 0 12 8 0 0 0 0 0 0 0 0 0 0 0 0 2 0 2 0 0 0 0 0 0 0 10 0 0 0 8 12 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 7 0 0 0 9 0 5 11 0 0 0 0 0 0 11 6 0 0 0 9 0 7 0 0 0 0 0 0 0 0 0 0 3 0 0 0 3 0 0 0 0 0 0 0 5 0 12 6 0 0 0 10 0 0 0 0 0 0 0 0 5 0 4 0 0 0 8 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]

"""
    findFacet(data, i, j, k, threshold, startIndex)

Return the position and index data of a mesh approximating the isosurface of `data` through the cell
in i, j, k coordinates (integer positions corresponding to the corners of grid cells)

"""
function findFacet(data::Voxel, i::Int, j::Int, k::Int, threshold, startIndex)

  corners = [
  i i i+1 i+1 i i i+1 i+1
  j j+1 j+1 j j j+1 j+1 j
  k k k k k+1 k+1 k+1 k+1
  ]

  data0 = data[corners[:,1]...]
  data1 = data[corners[:,2]...]
  data2 = data[corners[:,3]...]
  data3 = data[corners[:,4]...]
  data4 = data[corners[:,5]...]
  data5 = data[corners[:,6]...]
  data6 = data[corners[:,7]...]
  data7 = data[corners[:,8]...]
  cornerData = [data0, data1, data2, data3, data4, data5, data6, data7]

  #The vertex bitmask `index` is made by checking which corners lie above the threshold
  currVert = 0x01
  index = 0x00 :: UInt8
  for dat in cornerData
    if dat > threshold
      index |= currVert
    end
    currVert <<= 1
  end

  edges = edgeTable[index+1]

  vindex = startIndex
  edgeIndices = zeros(Int, 12)
  positions = nil(Vector{Float64})
  currEdge = 0x001
  #create the positions and indices for the positions of the isosurface on each intersected edge
  for j=1:12
    if edges & currEdge != 0
      edgeIndices[j] = vindex
      vertIndices = edgeCouplings[:,j]
      pos0 = corners[:,vertIndices[1]+1]
      pos1 = corners[:,vertIndices[2]+1]
      data0 = cornerData[vertIndices[1]+1]
      data1 = cornerData[vertIndices[2]+1]
      point = pos0 + (threshold - data0) * (pos1 - pos0) / (data1 - data0)
      positions = cons(point, positions)
      vindex += 1
    end
    currEdge <<= 1
  end

  polyVerts = triTable[:, index+1]
  indices = nil(Int)

  #using the generated positions and indices, write final sequence of indices
  #to the index array, according to the triangles making up the isosurface
  triCounter = 1
  for i=1:length(polyVerts)
    if polyVerts[i] == 0 break end
    if triCounter % 3 == 0
      indices = cons(edgeIndices[polyVerts[i-2]], indices)
      indices = cons(edgeIndices[polyVerts[i-1]], indices)
      indices = cons(edgeIndices[polyVerts[i]], indices)
    end
    triCounter += 1
  end
  return positions, indices
end

"concatenate the reverse of l1 to l2"
function revcat{T}(l1::LinkedList{T}, l2::LinkedList{T})
  for h in l1
    l2 = cons(h, l2)
  end
  return l2
end

function createMesh(grid::MeshGrid, threshold, startIndex=1, uniformScale=0.001)
  return createMesh(grid.data, (grid.maxPt - grid.minPt)..., grid.minPt..., threshold, startIndex, uniformScale)
end

"""
    createMesh(data, scaleX, scaleY, scaleZ, threshold, startIndex=1, uniformScale=3)

Create a mesh approximating the isosurface at level `threshold` of the scalar field `data`.
The scaling of the mesh is determined by `scaleX`, `scaleY`, and `scaleZ`
with a global scaling of `uniformScale`.
"""
function createMesh(data::Voxel, scaleX, scaleY, scaleZ, posX, posY, posZ, threshold, startIndex=1, uniformScale=0.001)
  positions = nil(Vector{Float64})
  indices = nil(Int)
  currIndex = startIndex
  dims = size(data)
  scaleX, scaleY, scaleZ = scaleX / dims[1] * uniformScale, scaleY / dims[2] * uniformScale, scaleZ / dims[3] * uniformScale
  posX, posY, posZ = posX * uniformScale, posY * uniformScale, posZ * uniformScale
  #compile the isosurfaces of all the cells in the volume into a single mesh
  for i=1:dims[1]-1
    for j=1:dims[2]-1
      for k=1:dims[3]-1
        pos, ind = findFacet(data, i, j, k, threshold, currIndex)
        positions = revcat(reverse(pos), positions)
        indices = revcat(reverse(ind), indices)
        currIndex += length(pos)
      end
    end
  end
  lenp = length(positions)
  leni = length(indices)
  positionsvec = Array(Float64, 3, lenp)
  indicesvec = Array(Int, leni)
  for (i, pos) in enumerate(positions)
    positionsvec[:,lenp-i+1] = pos .* [scaleX, scaleY, scaleZ] + [posX, posY, posZ]
  end
  for (i, ind) in enumerate(indices)
    indicesvec[leni-i+1] = ind
  end
  return positionsvec, indicesvec
end

type BoundsChecker
  disp::Vector{Float64}
  count::Int
  BoundsChecker() = new(zeros(3), 0)
end

type DoubleFinder
  index::Int
  indexMap::Dict{Int, Int}
  positions::Matrix{Float64}
  DoubleFinder() = new(1, Dict(), Array(Float64, 3, 0))
end

import Base.call
function call(finder::DoubleFinder, pos, ind)
  finder.positions = hcat(finder.positions, pos[:, 1])
  for i in ind
    finder.indexMap[i] = finder.index
  end
  finder.index += 1
end

function call(checker::BoundsChecker, pos, ind)
  disp = abs(maximum(pos, 2) - minimum(pos, 2))
  if disp[1] > checker.disp[1] checker.disp[1] = disp[1] end
  if disp[2] > checker.disp[2] checker.disp[2] = disp[2] end
  if disp[3] > checker.disp[3] checker.disp[3] = disp[3] end
  if size(pos)[2] > checker.count checker.count = size(pos)[2] end
end

"Return a list of int sets corresponding to the indices that now refer to the same position"
function findDoubles(positions, tol, verbose)
  finder = DoubleFinder()
  if verbose println("building octree...") end
  octree = Node(positions, minimum(positions, 2)[:,1] - [tol, tol, tol], maximum(positions, 2)[:,1] + [tol, tol, tol], 2, tol, [1:size(positions)[2];])
  if verbose
    println("points in octree: $(length(octree))")
    checker = BoundsChecker()
    map(checker, octree)
    println("maximum span of points in cells: $(checker.disp)")
    println("max number of points in cells: $(checker.count)")
    println("collecting doubles")
  end
  map(finder, octree)
  return finder
end

function removeDoubles(verts, indices, tol=1e-8, verbose=false)
  doubles = findDoubles(verts, tol, verbose)
  if verbose println("generating new mesh") end
  newIndices = Array(Int, 0)
  if verbose
    degenerates = 0
  end
  for j=1:div(length(indices), 3)
    index0 = indices[j*3 - 2]
    index1 = indices[j*3 - 1]
    index2 = indices[j*3]

    index0 = doubles.indexMap[index0]
    index1 = doubles.indexMap[index1]
    index2 = doubles.indexMap[index2]


    if index0 != index1 != index2 != index0
      newIndices = vcat(newIndices, [index0, index1, index2])
    elseif verbose
      degenerates += 1
    end
  end
  if verbose
    println("removed $(degenerates) degenerate faces")
    println("points in new mesh: $(size(doubles.positions)[2])")
  end
  return doubles.positions, newIndices
end

function countParts(maxIndex, indices)
  indexSet = IntSet(1:maxIndex)
  parts = 0
  while !isempty(indexSet)
    curr = first(indexSet)
    stack = Stack(Int)
    push!(stack, curr)
    delete!(indexSet, curr)
    found = false
    while !isempty(stack)
      curr = pop!(stack)
      for i=1:3:length(indices)
        if indices[i] == curr || indices[i+1] == curr || indices[i+2] == curr
          if in(indices[i], indexSet)
            delete!(indexSet, indices[i])
            push!(stack, indices[i])
            found = true
          end
          if in(indices[i+1], indexSet)
            delete!(indexSet, indices[i+1])
            push!(stack, indices[i+1])
            found = true
          end
          if in(indices[i+2], indexSet)
            delete!(indexSet, indices[i+2])
            push!(stack, indices[i+2])
            found = true
          end
        end
      end
    end
    if found parts += 1 end
  end
  return parts
end

function separate(indices)
  faceSet = IntSet(1:div(length(indices), 3))
  parts = nil(Vector{Int64})
  while !isempty(faceSet)
    curr = first(faceSet)
    stack = Stack(Int)
    push!(stack, curr)
    delete!(faceSet, curr)
    part = nil(Int)
    found = false
    while !isempty(stack)
      curr = pop!(stack)
      i0 = curr*3 - 2
      part = cons(indices[i0], part)
      part = cons(indices[i0+1], part)
      part = cons(indices[i0+2], part)
      for i=1:3:length(indices)
        f = div(i-1, 3)+1
        if in(f, faceSet)
          if indices[i] == indices[i0] || indices[i+1] == indices[i0] || indices[i+2] == indices[i0] || indices[i] == indices[i0+1] || indices[i+1] == indices[i0+1] || indices[i+2] == indices[i0+1] || indices[i] == indices[i0+2] || indices[i+1] == indices[i0+2] || indices[i+2] == indices[i0+2]
            push!(stack, f)
            delete!(faceSet, f)
            found = true
          end
        end
      end
    end
    if found
      len = length(part)
      partarr = Array(Int, len)
      for (i, ind) in enumerate(part)
        partarr[len-i+1] = ind
      end
      parts = cons(partarr, parts)
    end
  end
  return parts
end
