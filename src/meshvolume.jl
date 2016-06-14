signedvolume(p1, p2, p3) = dot(cross(p1, p2), p3)/6.0

"Compute the volume of a mesh defined by `points` and `indices`. Mesh must be closed with outward facing normals."
function volume(points, indices)
  vol = 0
  for i=1:div(length(indices), 3)
    p1 = points[:, indices[(i-1)*3 + 1]]
    p2 = points[:, indices[(i-1)*3 + 2]]
    p3 = points[:, indices[(i-1)*3 + 3]]
    vol += signedvolume(p1, p2, p3)
  end
  return vol
end
