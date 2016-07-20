const xplus = 1
const xminus = -1
const yplus = 2
const yminus = -2
const zplus = 3
const zminus = -3
const mindisp = 0.1

lerp(x, y, alpha) = x*alpha + y*(1-alpha)

function linearFilter(filter, a::Vector)
    len = length(a)
    b = zeros(len)
    for i=1:len
        for j=1:len
            b[i] += filter(j - i) * a[j]
        end
    end
    return b
end

function makeGaussianFilter(stdev)
    g(x) = exp(-x^2/(2 * stdev^2))
    f(a) = linearFilter(g, a)
    return f
end

function postoarray(dims, minPt::Vector, maxPt::Vector, pos::Vector)
    disp = maxPt - minPt
    relPos = pos - minPt
    realPos = relPos ./ disp .* [dims...]
    realx = clamp(ceil(Int, realPos[1]), 1, dims[1])
    realy = clamp(ceil(Int, realPos[2]), 1, dims[2])
    realz = clamp(ceil(Int, realPos[3]), 1, dims[3])
    realInd = [realx, realy, realz]
    offset = realPos - realInd + 0.5
    return realInd, offset
end

"A uniform grid containing density information loaded from point data"
immutable MeshGrid
    minPt::Array{Float64, 1}
    maxPt::Array{Float64, 1}
    data::Array{Float64, 3}

    """
    MeshGrid(x_resolution, y_resolution, z_resolution, hits, dilate=false, stdev=20.0, buffer=[0, 0, 0], threshold=Inf)

    Read [x, y, z, E] from the columns of `hits` and build a 3D density volume.
    If `dilate` is enabled, apply a gaussian filter with radius `stdev`.
    `buffer` controls the expansion of the region beyond the bounding box of the hits (where 1 is double the dimensions.)
    """
    MeshGrid(resx, resy, resz, hits, dilate=false, stdev=20.0, buffer = [0.0, 0.0, 0.0], threshold=Inf) = MeshGrid(minimum(hits, 2)[1:3], maximum(hits, 2)[1:3], resx, resy, resz, hits, dilate, stdev, buffer, threshold)

    function MeshGrid(mn::Vector, mx::Vector, resx::Int, resy::Int, resz::Int, hits::Matrix, dilate::Bool, stdev::Real, buffer::Vector, threshold::Real)
        if stdev == 0 dilate = false end
        initdisp = mx - mn
        for i=1:3
            if initdisp[i] == 0
                initdisp[i] = mindisp
                mx[i] += mindisp
            end
        end
        recompute = true
        minPt = nothing
        maxPt = nothing
        data = nothing
        while recompute
            minPt = mn - initdisp .* buffer / 2
            maxPt = mx + initdisp .* buffer / 2
            disp = maxPt - minPt
            dims = [resx, resy, resz]
            data = zeros(resx, resy, resz)

            #first pass: read the position data directly into the corresponding cells
            for i=1:size(hits, 2)
                pos = hits[1:3, i]
                data[postoarray(dims, minPt, maxPt, pos)[1]...] += hits[4, i]

            end
            #second pass: blur the data (optional)
            if dilate
                factor = dims ./ disp
                stdevx, stdevy, stdevz = ((stdev * factor)...)

                gaussianFilters = (makeGaussianFilter(stdevx), makeGaussianFilter(stdevy), makeGaussianFilter(stdevz))
                for j=1:3
                    data = mapslices(gaussianFilters[j], data, j)
                end
            end

            #recompute
            if threshold < Inf
                maxEnergyDensity = maximum(data)
                minEnergyDensity = minimum(data)
                maxima = [maximum(data[1,:,:]), maximum(data[end,:,:]),
                maximum(data[:,1,:]), maximum(data[:,end,:]),
                maximum(data[:,:,1]), maximum(data[:,:,end])]
                boundmaxx = maximum(maxima[1:2])
                boundmaxy = maximum(maxima[3:4])
                boundmaxz = maximum(maxima[5:6])
                if boundmaxx >= threshold
                    buffer[1] *= 2
                end
                if boundmaxy >= threshold
                    buffer[2] *= 2
                end
                if boundmaxz >= threshold
                    buffer[3] *= 2
                end
                if maximum([boundmaxx, boundmaxy, boundmaxz]) < threshold
                    recompute = false
                    print("border: $(buffer)")
                else
                    print("(resizing)")
                end

            else
                recompute = false
            end
        end
        new(minPt, maxPt, data)
    end
end

import Base.getindex
"Get the energy of the field at position x, y, z"
function getindex(grid::MeshGrid, x::Real, y::Real, z::Real)
    dims = size(grid.data)
    realInd, offset = postoarray(dims, grid.minPt, grid.maxPt, [x, y, z])

    xother = clamp(offset[1] > 0 ? realInd[1] + 1 : realInd[1] - 1, 1, dims[1])
    yother = clamp(offset[2] > 0 ? realInd[2] + 1 : realInd[2] - 1, 1, dims[2])
    zother = clamp(offset[3] > 0 ? realInd[3] + 1 : realInd[3] - 1, 1, dims[3])

    offset = abs(offset)

    v000 = grid.data[realInd...]
    v100 = grid.data[xother,     realInd[2], realInd[3]]
    v010 = grid.data[realInd[1], yother,     realInd[3]]
    v110 = grid.data[xother,     yother,     realInd[3]]
    v001 = grid.data[realInd[1], realInd[2], zother]
    v101 = grid.data[xother,     realInd[2], zother]
    v011 = grid.data[realInd[1], yother,     zother]
    v111 = grid.data[xother,     yother,     zother]

    v00 = lerp(v001, v000, offset[3])
    v10 = lerp(v101, v100, offset[3])
    v01 = lerp(v011, v010, offset[3])
    v11 = lerp(v111, v110, offset[3])

    v0 = lerp(v01, v00, offset[2])
    v1 = lerp(v11, v10, offset[2])

    v = lerp(v1, v0, offset[1])
    return v
end

# "Get the energy of the cell at `i`, `j`, `k` using the density and volume of the cell."
# function getenergy(grid::MeshGrid, i, j, k)
#   dims = size(grid.data)
#   cellvol = prod((grid.maxPt - grid.minPt) ./ [dims[1], dims[2], dims[3]]) :: Float64
#   return grid.data[i, j, k] * cellvol
# end

"Compute the edge weight from the given voxel along the given axis (as defined in this module) using the provided binary function func of the pixel values"
function edgeWeight(grid::MeshGrid, func, i, j, k, axis)
    if axis == xplus
        other = grid.data[i+1, j, k]
    elseif axis == xminus
        other = grid.data[i-1, j, k]
    elseif axis == yplus
        other = grid.data[i, j+1, k]
    elseif axis == yminus
        other = grid.data[i, j-1, k]
    elseif axis == zplus
        other = grid.data[i, j, k+1]
    elseif axis == zminus
        other = grid.data[i, j, k-1]
    else
        println("unexpected axis $(axis) for edge weight")
        return 0
    end

    return func(grid.data[i, j, k], other)
end

immutable Edge
    x::Int
    y::Int
    z::Int
    dir::Int
    weight::Float64
    #Edge(x, y, z, dir) = new(x, y, z, dir)
end

"get the point at the end of the edge"
function getendpoint(e::Edge)
    if e.dir == xplus
        return (e.x+1, e.y, e.z)
    elseif e.dir == yplus
        return (e.x, e.y+1, e.z)
    elseif e.dir == zplus
        return (e.x, e.y, e.z+1)
    else
        return (-1, -1, -1)
    end
end

"""
segment(grid, k_int=1.0, func=(a, b) -> abs(a-b))

Partition the set of all grid coordinates into disjoint subsets, satisfying the courseness and fineness
criterea specified in the penalty function
"""
function segment(grid::MeshGrid, kint=1.0, func=(a, b)->abs(a-b))
    edges = Array(Edge, 0)
    dims = size(grid.data)
    subsets = Array(Tuple{Float64, Array{Tuple{Int, Int, Int}, 1}}, 0) #an array of disjoint subsets, represented by an array of tuples of (integer) positions and an internal difference value.

    for i=1:dims[1]
        for j=1:dims[2]
            for k=1:dims[3]
                #if i < dims[1] && j < dims[2] && k < dims[3]
                if i < dims[1]
                    weightx = edgeWeight(grid, func, i, j, k, xplus)
                    push!(edges, Edge(i, j, k, xplus, weightx))
                end
                if j < dims[2]
                    weighty = edgeWeight(grid, func, i, j, k, yplus)
                    push!(edges, Edge(i, j, k, yplus, weighty))
                end
                if k < dims[3]
                    weightz = edgeWeight(grid, func, i, j, k, zplus)
                    push!(edges, Edge(i, j, k, zplus, weightz))
                end

                #end
                push!(subsets, (0.0, [(i, j, k)]))
            end
        end
    end
    sort(edges, by=edge->edge.weight)
    for (id, edge) in enumerate(edges)
        #find subsets in which the neighboring nodes currently reside
        indices1 = find(map(ss->in((edge.x, edge.y, edge.z), ss[2]), subsets))
        indices2 = find(map(ss->in(getendpoint(edge), ss[2]), subsets))

        if length(indices1) == 0
            continue
        end
        if length(indices2) == 0
            continue
        end

        index1 = indices1[1]
        index2 = indices2[1]
        if index1 == index2
            continue #the edge is not in disjoint subsets
        end
        subset1 = subsets[index1]
        subset2 = subsets[index2]

        #compute their mutual min internal difference
        mint = min(subset1[1] + kint/length(subset1[2]), subset2[1] + kint/length(subset2[2]))

        if edge.weight < mint
            #println("merging subset $(index1) and $(index2) with mInt = $(mint) and k=$(k)")
            subsets[index1] = (edge.weight, vcat(subset1[2], subset2[2]))
            deleteat!(subsets, index2)
        end
    end
    return subsets
end
