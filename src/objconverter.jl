using Formatting
"saveObj(positions, indices, \"filename\") saves the mesh represented by the given verts and tris to a file"
function saveObj(positions, indices, filename="untitled.obj")
  f = open(filename, "w")
  for i=1:size(positions)[2]
    pos = positions[:,i]
    write(f, "v $(pos[1]) $(pos[2]) $(pos[3])\n")
  end
  for i=1:div(length(indices), 3)
    write(f, "f $(indices[i*3 - 2]) $(indices[i*3 - 1]) $(indices[i*3])\n")
  end
  close(f)
end

function saveOff(positions, indices, filename="untitled.obj")
  f = open(filename, "w")
  write(f, "OFF\n")
  write(f, "$(length(positions)) $(div(length(indices), 3)) $(edgecount(indices))\n")
  for i=1:size(positions)[2]
    pos = positions[:,i]
    printfmt(f, "{1:.6f} {2:.6f} {3:.6f}\n", pos[1], pos[2], pos[3])
  end
  for i=1:div(length(indices), 3)
    write(f, "3  $(indices[i*3 - 2]-1) $(indices[i*3 - 1]-1) $(indices[i*3]-1)\n")
  end
  close(f)
end

function edgecount(indices)
    edgeSet = Set{Vector{Int}}()
    for i=1:div(length(indices), 3)
        push!(edgeSet, sort([indices[i*3-2], indices[i*3-1]]))
        push!(edgeSet, sort([indices[i*3-1], indices[i]]))
        push!(edgeSet, sort([indices[i], indices[i*3-2]]))
    end
    return length(edgeSet)
end
