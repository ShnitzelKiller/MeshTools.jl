using DataStructures
#include("objconverter.jl")
#f0 is the face bordering the edge between i0 and i1
normalize(x) = x/norm(x)
type Face
  i0::Int
  i1::Int
  i2::Int
  conflicts::LinkedList{Int}
  disabled::Bool
  farthest::Tuple{Int, Float64}
  f0::Face
  f1::Face
  f2::Face
  area::Float64
  normal::Vector{Float64}
  Face(i0, i1, i2) = new(i0, i1, i2, nil(Int), false, (-1, -1.0))
end

import Base.show
show(f::Face) = "($(f.i0), $(f.i1), $(f.i2))"

function computeface(face::Face, points::Matrix)
  p0 = points[:,face.i0]
  p1 = points[:,face.i1]
  p2 = points[:,face.i2]
  d = cross(p1-p0, p2-p0)
  area = length(d)
  face.normal = d/area
  face.area = area/2.0
end

type Mesh
  faces::LinkedList{Face}
  vertices::Matrix
  Mesh(points) = new(nil(Face), points)
end

import Base.push!
function push!(mesh::Mesh, face::Face)
  mesh.faces = cons(face, mesh.faces)
end

function getmeshdata(mesh::Mesh, simplify=true)
  indices = nil(Int)
  for face in mesh.faces
    if !face.disabled
      indices = cons(face.i0, indices)
      indices = cons(face.i1, indices)
      indices = cons(face.i2, indices)
    end
  end
  indices = reverse(indices)
  if simplify
    used = falses(size(mesh.vertices)[2])
    for i in indices
      used[i] = true
    end
    usedlength = 0
    for u in used
      if u usedlength += 1 end
    end
    indexmap = zeros(size(mesh.vertices)[2])
    curr = 1
    for i=1:length(indexmap)
      indexmap[i] = curr
      if used[i]
        curr += 1
      end
    end
    newindices = Array(Int, length(indices))
    for (i, ind) in enumerate(indices)
      newindices[i] = indexmap[ind]
    end
    newpositions = Array(Float64, 3, usedlength)
    curr = 1
    for i=1:size(mesh.vertices)[2]
      if used[i]
        newpositions[:,curr] = mesh.vertices[:,i]
        curr += 1
      end
    end
    return newpositions, newindices
  else
    indexarray = Array(Int, length(indices))
    for (i, ind) in enumerate(indices)
      indexarray[i] = ind
    end
    return mesh.vertices, indexarray
  end
end

function computeconflictlists(mesh::Mesh, i::Int)
  maxDist = 0
  farFace = head(mesh.faces)
  for face in mesh.faces
    if !face.disabled
      p0 = mesh.vertices[:,face.i0]
      pos = mesh.vertices[:,i]
      dist = dot(pos-p0, face.normal)
      if dist > maxDist
        maxDist = dist
        farFace = face
      end
    end
  end
  if maxDist > 0
    farFace.conflicts = cons(i, farFace.conflicts)
    if maxDist > farFace.farthest[2]
      farFace.farthest = (i, maxDist)
    end
  end
end

function getfarthestpoints(points::Matrix, debug=false)
  n = size(points)[2]

  i0, i1, i2, i3 = -1, -1, -1, -1
  maxDist = 0.0
  for i=1:n
    for j=i+1:n
      dist = norm(points[:,i] - points[:,j])
      if dist > maxDist
        maxDist = dist
        i0 = i
        i1 = j
      end
    end
  end

  if i0 == -1
    return 1, 2, 3, 4
  elseif debug
    println("first two points: $(points[:,i0]) at $i0 and $(points[:,i1]) at $i1")
  end
  if i1 < i0 i1, i0 = i0, i1 end
  disp = points[:,i1] - points[:,i0]
  maxDist = 0.0
  for i=[1:i0-1; i0+1:i1-1; i1+1:n]
    currdisp = points[:,i] - points[:,i0]
    proj = disp * dot(disp, currdisp) / dot(disp, disp)
    dist = norm(currdisp - proj)
    if dist > maxDist
      i2 = i
      maxDist = dist
    end
  end

  #if any index wasn't found, just default to the safe 1, 2, 3, 4
  if i2 == -1
    return 1, 2, 3, 4
  elseif debug
    println("third point: $(points[:,i2]) at $i2")
  end



  i0, i1, i2 = sort([i0, i1, i2])
  disp = points[:,i1] - points[:,i0]
  disp2 = points[:,i2] - points[:,i0]
  maxDist = 0.0
  for i=[1:i0-1; i0+1:i1-1; i1+1:i2-1; i2+1:n]
    normal = normalize(cross(disp, disp2))
    currdisp = points[:,i] - points[:,i0]
    dist = abs(dot(currdisp, normal))
    if dist > maxDist
      i3 = i
      maxDist = dist
    end
  end

  #if any index wasn't found, just default to the safe 1, 2, 3, 4
  if i3 == -1
    return 1, 2, 3, 4
  elseif debug
    println("fourth point: $(points[:,i3]) at $i3")
  end

  if dot(points[:,i3] - points[:,i0], cross(points[:,i1] - points[:,i0], points[:,i2] - points[:,i0])) > 0
    i1, i2 = i2, i1
  end
  return i0, i1, i2, i3
end

"""
    convexhull(points, simplify=true, epsilon=1e-8, debug=false)

construct the convex hull of this mesh using the 3D points specified in the
columns of `points`.

"""
function convexhull(points::Matrix, simplify=true, epsilon=0, debug=false)



  #first, build the first tetrahedron using the widest pair of points, and the farthest point from the line, then from the plane
  n = size(points)[2]
  if n < 4
    println("too few points (minimum 4)")
    return [], []
  end

  i0, i1, i2, i3 = getfarthestpoints(points, debug)

  face012 = Face(i0, i1, i2)
  face023 = Face(i0, i2, i3)
  face321 = Face(i3, i2, i1)
  face031 = Face(i0, i3, i1)

  computeface(face012, points)
  computeface(face023, points)
  computeface(face321, points)
  computeface(face031, points)

  face012.f0 = face031
  face012.f1 = face321
  face012.f2 = face023

  face023.f0 = face012
  face023.f1 = face321
  face023.f2 = face031

  face321.f0 = face023
  face321.f1 = face012
  face321.f2 = face031

  face031.f0 = face023
  face031.f1 = face321
  face031.f2 = face012


  mesh = Mesh(points)
  push!(mesh, face012)
  push!(mesh, face023)
  push!(mesh, face321)
  push!(mesh, face031)

  i0, i1, i2, i3 = sort([i0, i1, i2, i3])
  if debug
    println("start: $i0, $i1, $i2, $i3")
  end
  #partition vertex indices into conflict lists
  for i=[1:i0-1; i0+1:i1-1; i1+1:i2-1; i2+1:i3-1; i3+1:n]
    computeconflictlists(mesh, i)
  end

  step = 0
  while true
    if debug
      assertconflicts(mesh, step)
      assertmesh(mesh)
      step += 1
    end

    farthest, maxDist = -1, 0.0
    farFace = head(mesh.faces)
    for face in mesh.faces
      if (!face.disabled && face.farthest[2] > maxDist)
        farthest, maxDist = face.farthest
        farFace = face
      end
    end

    if farthest != -1
      frontier = Stack(Face)
      visited = nil(Face)
      push!(frontier, farFace)
      farFace.disabled = true
      horizon = nil(Tuple{Int, Int})
      outsidefaces = nil(Face)

      #BFS of visible faces
      while !isempty(frontier)
        face = pop!(frontier)
        visited = cons(face, visited)
        for (i, nextface) in [(1, face.f0), (2, face.f1), (3, face.f2)]
          #compute index into previous face so we start at the next index in CCW order
          if nextface.disabled continue end

          #visibility check
          v = points[:,farthest] - points[:,nextface.i0]
          if dot(v, nextface.normal) > epsilon
            push!(frontier, nextface)
            nextface.disabled = true
          else
            #add to the list of horizon edges
            # println("non-visible face: $(show(nextface))")
            edges = [face.i0, face.i1, face.i2]
            horizon = cons((edges[i], edges[(i%3) + 1]), horizon)
            outsidefaces = cons(nextface, outsidefaces)
          end
        end
      end

      #create new faces along horizon
      len = length(horizon)
      # println("length of horizon: $(len)")
      # println("horizon: $(horizon)")
      # println("Farthest: $(farthest)")
      # println("number of faces visited: $(length(visited))")
      horizon, outsidefaces = makesequential(horizon, outsidefaces)
      # lastpoint = first(horizon)[2]
      # cyclic=true
      # for edge in horizon
      #   if lastpoint != edge[2]
      #     cyclic=false
      #     break
      #   end
      #   lastpoint = edge[1]
      # end
      # if !cyclic
      #   println("non-cyclic horizon: $(horizon), adding point $(points[:,farthest])")
      # else
      #   println("cyclic horizon: $(horizon)")
      # end
      newfaces = Array(Face, len)
      for (i, edge, outface) in zip(1:len, horizon, outsidefaces)
        newface = Face(edge[1], edge[2], farthest)
        computeface(newface, points)
        newface.f0 = outface
        newfaces[i] = newface
        push!(mesh, newface)
        if outface.i0 == newface.i1
          outface.f0 = newface
        elseif outface.i1 == newface.i1
          outface.f1 = newface
        elseif outface.i2 == newface.i1
          outface.f2 = newface
        else
          @assert false "outer face $(show(outface)) not bound (may cause holes)"
        end
      end
      #properly assign new neighbors
      for i=1:len
        j = i-1
        if j < 1
          j = len
        end
        k = i+1
        if k > len
          k = 1
        end
        newfaces[i].f1 = newfaces[j]
        newfaces[i].f2 = newfaces[k]
      end

      #distribute vertices (except the one being added to the hull) among conflict lists of enabled faces
      if debug assertunique(visited) end
      for face in visited
        for ind in face.conflicts
          if ind != farthest
            computeconflictlists(mesh, ind)
          end
        end
      end
    else
      return getmeshdata(mesh, simplify)
    end
  end
end

function assertunique{T}(things::LinkedList{T})
  thingset = Dict{T, T}()
  for thing in things
    @assert !haskey(thingset, thing) "duplicate $T found in list, $(thingset[thing]===thing ? "identical" : "not indentical")"
    thingset[thing] = thing
  end
end

function makesequential(edges, outfaces)
  len = length(edges)
  seq = list(first(edges))
  seq2 = list(first(outfaces))
  for i=1:len-1
    found = false
    for (edge, face) in zip(edges, outfaces)
      if edge[1] == first(seq)[2]
        seq = cons(edge, seq)
        seq2 = cons(face, seq2)
        found = true
        break
      end
    end
    @assert found "non-sequential edges in horizon: $(edges)"
  end
  return seq, seq2
end

function assertmesh(mesh::Mesh)
  for face in mesh.faces
    if face.disabled continue end
    neighb = [face.f0, face.f1, face.f2]
    inds = [face.i0, face.i1, face.i2]
    for i=1:3
      if neighb[i].f0 === face
        startind = 1
      elseif neighb[i].f1 === face
        startind = 2
      elseif neighb[i].f2 === face
        startind = 3
      else
        startind = 0
      end
      @assert startind != 0 "asymmetric conjoinment: $(show(neighb[i])) and $(show(face))"
      errstr = "incorrect conjoinment: $(show(neighb[i])) and $(show(face)) conjoined at index $startind and $i, $(neighb[i]===face ? "faces identical" : "faces different")"
      if startind == 1
        @assert (neighb[i].i0 == inds[i%3 + 1] && neighb[i].i1 == inds[i]) errstr
      elseif startind == 2
        @assert (neighb[i].i1 == inds[i%3 + 1] && neighb[i].i2 == inds[i]) errstr

      elseif startind == 3
        @assert (neighb[i].i2 == inds[i%3 + 1] && neighb[i].i0 == inds[i]) errstr
      end
    end
  end
end

function assertconflicts(mesh::Mesh, step::Int)
  indset = Set{Int}()
  dups = 0
  dupdict = Dict{Int, Int}()
  duplist = nil(Int)
  for f in mesh.faces
    if f.disabled continue end
    for ind in f.conflicts
      if in(ind, indset)
        dups += 1
        dupdict[ind] += 1
        duplist = cons(ind, duplist)
      else
        dupdict[ind] = 1
        push!(indset, ind)
      end
    end
  end
  dict2 = Dict{Int, Int}()
  for dup in duplist
    dict2[dup] = dupdict[dup]
  end
  @assert dups == 0 "$dups duplicates at step $step: $(dict2)"
end
