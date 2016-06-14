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
