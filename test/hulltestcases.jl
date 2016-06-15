using MeshTools

using Base.Test

function assertvalid(v, ind)
  len = size(v)[2]
  for i in ind
    @assert i <= len || i <= 1 "invalid mesh index"
  end
end

#test adding inside point: if error, it's seeing all faces (inside out)
data = hcat([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0.2, 0.2, 0.2])
v, ind = convexhull(data, true, 1e-8, true)
assertvalid(v, ind)
@test size(v)[2] == 4
@test length(ind) == 12

#adding an extra vertex
data = hcat([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1])
v, ind = convexhull(data, true, 1e-8, true)
@test size(v)[2] == 5
@test length(ind) == 18

#occluding a former hull vertex
data = hcat([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [2, 0, 0])
v, ind = convexhull(data, true, 1e-8, true)
@test size(v)[2] == 4
@test length(ind) == 12

v, ind = convexhull([1 2 3 4; 0 0 0 0; 0 0 0 0])
@test size(v)[2] == 4
@test length(ind) == 12

for i=1:20
  data = randn(3, 1000)
  v, ind = convexhull(data, true, 1e-8, true)
  assertvalid(v, ind)
end
