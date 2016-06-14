include("hull.jl")
include("objconverter.jl")
include("meshvolume.jl")
#using PyPlot

data = hcat([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1], [0, 2, 0.5], [1, 2, 0.5])
#data = hcat(data, [1, 1.1, 2.2])

# data = randn(3, 500)
# data = hcat(extra, data)

#PyPlot.scatter3D(data[1,:], data[2,:], data[3,:])
#plt[:show]()

v, ind = convexhull(data, true)
saveObj(v, ind, "../output/rand.obj")
println("volume of hull: $(volume(v, ind))")
