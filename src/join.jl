function joinmeshes(meshes...)
    len, indlen = foldl((x, y)->(x[1]+size(y[1], 2), x[2]+length(y[2])), (0, 0), meshes)
    vertices = Array(Float64, 3, len)
    indices = Array(Int, indlen)
    offset = 0
    indoffset = 0
    for (v, ind) in meshes
        l = size(v, 2)
        for i=1:l
            vertices[:, offset+i] = v[:, i]
        end
        for (i, index) in enumerate(ind)
            indices[indoffset+i] = index + offset
        end
        offset += l
        indoffset += length(ind)
    end
    return vertices, indices
end
