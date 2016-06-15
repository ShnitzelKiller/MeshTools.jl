using MeshTools

data1 = randn(3, 20)
data2 = randn(3, 20)
data = hcat(data1, [100, 0, 0], data2)

v, ind = convexhull(data, true, 1e-8, true)
println("volume of hull: $(volume(v, ind))")
saveObj(v, ind, "output/rand.obj")
