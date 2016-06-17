signedvolume(p1, p2, p3) = dot(cross(p1, p2), p3)/6.0

"Compute the volume of a mesh defined by `points` and `indices`. Mesh must be closed with outward facing normals."
function volume(points, indices)
  vol = 0
  meanpoint = mean(points, 2)[:, 1]
  for i=1:div(length(indices), 3)
    p1 = points[:, indices[(i-1)*3 + 1]] - meanpoint
    p2 = points[:, indices[(i-1)*3 + 2]] - meanpoint
    p3 = points[:, indices[(i-1)*3 + 3]] - meanpoint
    vol += signedvolume(p1, p2, p3)
  end
  return vol
end

"Compute the surface area of a mesh defined by `points` and `indices`."
function surfacearea(points, indices)
  surf = 0
  for i=1:div(length(indices), 3)
    p1 = points[:, indices[(i-1)*3 + 1]]
    p2 = points[:, indices[(i-1)*3 + 2]]
    p3 = points[:, indices[(i-1)*3 + 3]]
    v2 = p2 - p1
    v3 = p3 - p1
    area = norm(cross(v2, v3))/2.0
    surf += area
  end
  return surf
end
