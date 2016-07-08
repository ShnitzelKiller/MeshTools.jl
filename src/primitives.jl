function cylinder(r, h, n)
    v = Array(Float64, 3, 2*n)
    ind = Array(Int, 6*n)
    for i=1:n
        p1 = [h/2, r*cos(i/n*2pi), r*sin(i/n*2pi)]
        p2 = [-h/2, r*cos(i/n*2pi), r*sin(i/n*2pi)]
        v[:, i*2-1] = p1
        v[:, i*2] = p2
        ind[i*6-5] = i*2-1
        ind[i*6-4] = i*2
        ind[i*6-3] = mod(i*2, n*2) + 1
        ind[i*6-2] = mod(i*2, n*2) + 1
        ind[i*6-1] = i*2
        ind[i*6]   = mod(i*2+1, n*2) + 1
    end
    return v, ind
end
