immutable Node
  data::Nullable{Matrix}
  info::Nullable{Vector}
  isLeaf::Bool

  minPt::Vector
  maxPt::Vector

  child000::Node
  child001::Node
  child010::Node
  child011::Node
  child100::Node
  child101::Node
  child110::Node
  child111::Node

"""
   Node(data, minPt, maxPt, splitThreshold, tol, info)

 Create an octree using the points in the colums of `data`, the specified min and max points,
 and maximum number of points per node `splitThreshold`, with some extra information associated with
 each point, `info`. The depth is limited by the minimum cell width `tol`.
"""
  function Node{T}(data::Matrix, minPt::Vector, maxPt::Vector, splitThreshold::Int, tol, info::Vector{T})

    if size(data)[2] == 0
      new(Nullable{Matrix}(), Nullable{Vector}(), true)
    elseif size(data)[2] > splitThreshold && maximum(abs(maxPt - minPt)) > tol * 2
      midPt = (minPt + maxPt)/2.0

      data000 = Array(Float64, 3, 0)
      data001 = Array(Float64, 3, 0)
      data010 = Array(Float64, 3, 0)
      data011 = Array(Float64, 3, 0)
      data100 = Array(Float64, 3, 0)
      data101 = Array(Float64, 3, 0)
      data110 = Array(Float64, 3, 0)
      data111 = Array(Float64, 3, 0)

      info000 = Array(T, 0)
      info001 = Array(T, 0)
      info010 = Array(T, 0)
      info011 = Array(T, 0)
      info100 = Array(T, 0)
      info101 = Array(T, 0)
      info110 = Array(T, 0)
      info111 = Array(T, 0)

      for i=1:size(data)[2]
        point = data[:,i]
        if point[1] < midPt[1]
          if point[2] < midPt[2]
            if point[3] < midPt[3]
              data000 = hcat(data000, point)
              info000 = vcat(info000, info[i])
            else
              data001 = hcat(data001, point)
              info001 = vcat(info001, info[i])
            end
          else
            if point[3] < midPt[3]
              data010 = hcat(data010, point)
              info010 = vcat(info010, info[i])
            else
              data011 = hcat(data011, point)
              info011 = vcat(info011, info[i])
            end
          end
        else
          if point[2] < midPt[2]
            if point[3] < midPt[3]
              data100 = hcat(data100, point)
              info100 = vcat(info100, info[i])
            else
              data101 = hcat(data101, point)
              info101 = vcat(info101, info[i])
            end
          else
            if point[3] < midPt[3]
              data110 = hcat(data110, point)
              info110 = vcat(info110, info[i])
            else
              data111 = hcat(data111, point)
              info111 = vcat(info111, info[i])
            end
          end
        end
      end
      child000 = Node(data000, minPt, midPt, splitThreshold, tol, info000)
      child001 = Node(data001, [minPt[1], minPt[2], midPt[3]], [midPt[1], midPt[2], maxPt[3]], splitThreshold, tol, info001)
      child010 = Node(data010, [minPt[1], midPt[2], minPt[3]], [midPt[1], maxPt[2], midPt[3]], splitThreshold, tol, info010)
      child011 = Node(data011, [minPt[1], midPt[2], midPt[3]], [midPt[1], maxPt[2], maxPt[3]], splitThreshold, tol, info011)
      child100 = Node(data100, [midPt[1], minPt[2], minPt[3]], [maxPt[1], midPt[2], midPt[3]], splitThreshold, tol, info100)
      child101 = Node(data101, [midPt[1], minPt[2], midPt[3]], [maxPt[1], midPt[2], maxPt[3]], splitThreshold, tol, info101)
      child110 = Node(data110, [midPt[1], midPt[2], minPt[3]], [maxPt[1], maxPt[2], midPt[3]], splitThreshold, tol, info110)
      child111 = Node(data111, midPt, maxPt, splitThreshold, tol, info111)

      new(Nullable{Matrix}(), Nullable{Vector}(), false, minPt, maxPt, child000, child001, child010, child011, child100, child101, child110, child111)
    else
      new(Nullable(data), Nullable(info), true)
    end
  end

end
import Base.map
"Apply the callable f to all non-null leaf nodes"
function map(f, octree::Node)
  if octree.isLeaf
    if !isnull(octree.data)
      f(get(octree.data), get(octree.info))
    end
  else
    map(f, octree.child000)
    map(f, octree.child001)
    map(f, octree.child010)
    map(f, octree.child011)
    map(f, octree.child100)
    map(f, octree.child101)
    map(f, octree.child110)
    map(f, octree.child111)
  end
end

type VertexCounter
  count::Int
  VertexCounter() = new(0)
end

@compat (function (counter::VertexCounter)(pos, ind)
  counter.count += length(ind)
end)

import Base.length
function length(octree::Node)
  counter = VertexCounter()
  map(counter, octree)
  return counter.count
end

import Base.in

function doAt(f, point::Vector, octree::Node, default)
  if octree.isLeaf
    if !isnull(octree.data)
      data = get(octree.data)
      info = get(octree.info)
      return f(data, info)
    else
      return default
    end
  else
    midPt = (octree.minPt + octree.maxPt)/2.0
    if point[1] < midPt[1]
      if point[2] < midPt[2]
        if point[3] < midPt[3]
          return doAt(f, point, octree.child000, default)
        else
          return doAt(f, point, octree.child001, default)
        end
      else
        if point[3] < midPt[3]
          return doAt(f, point, octree.child010, default)
        else
          return doAt(f, point, octree.child011, default)
        end
      end
    else
      if point[2] < midPt[2]
        if point[3] < midPt[3]
          return doAt(f, point, octree.child100, default)
        else
          return doAt(f, point, octree.child101, default)
        end
      else
        if point[3] < midPt[3]
          return doAt(f, point, octree.child110, default)
        else
          return doAt(f, point, octree.child111, default)
        end
      end
    end
  end
end

function in(point::Vector, inf, octree::Node)
  function f(data, info)
    for i=1:size(data)[2]
      if data[:,i] == point && info[i] == inf
        return true
      end
    end
    return false
  end
  return doAt(f, point, octree, false)
end

import Base.getindex
function getindex(octree::Node, x, y, z)
  point = [x, y, z]
  function f(data, info)
    infos = info[1:0]
    for i=1:size(data)[2]
      if data[:,i] == point
        push!(infos, info[i])
      end
    end
    return infos
  end
  return doAt(f, point, octree, [])
end
