include("hull.jl")

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
@assert size(v)[2] == 4 "$(size(v)[2]) vertices instead of 4"
@assert length(ind) == 4 * 3 "$(length(ind)) indices instead of 12 (3 for 4 faces)"

#adding an extra vertex
data = hcat([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1])
v, ind = convexhull(data, true, 1e-8, true)
@assert size(v)[2] == 5 "$(size(v)[2]) vertices instead of 5"
@assert length(ind) == 18 "$(length(ind)) indices instead of 18 (3 for 6 faces)"

#occluding a former hull vertex
data = hcat([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [2, 0, 0])
v, ind = convexhull(data, true, 1e-8, true)
@assert size(v)[2] == 5 "$(size(v)[2]) vertices instead of 5"
@assert length(ind) == 18 "$(length(ind)) indices instead of 18 (3 for 6 faces)"

for i=1:20
  data = randn(3, 50)
  v, ind = convexhull(data, true, 1e-8, true)
  assertvalid(v, ind)
end
