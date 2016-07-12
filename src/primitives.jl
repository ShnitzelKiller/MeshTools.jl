function cylinder(r, h, n)
    v = Array(Float64, 3, 2*n)
    ind = Array(Int, 6*n)
    for i=1:n
        p1 = [r*cos(i/n*2pi), h/2, r*sin(i/n*2pi)]
        p2 = [r*cos(i/n*2pi), -h/2, r*sin(i/n*2pi)]
        v[:, i*2-1] = p1
        v[:, i*2] = p2
        ind[i*6-5] = i*2-1
        ind[i*6-4] = mod(i*2, n*2) + 1
        ind[i*6-3] = i*2
        ind[i*6-2] = mod(i*2, n*2) + 1
        ind[i*6-1]   = mod(i*2+1, n*2) + 1
        ind[i*6] = i*2
    end
    return v, ind
end
