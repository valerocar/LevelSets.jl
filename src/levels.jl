#
# Level Sets
#
#    Functions for computing level curves and surfaces in simplicial complexes
#
#    Author: Carlos Valero Valdés
#    e-mail: valerocar@gmail.com
#
#

"""
    levelcurve(ps::Array{<:Real,2}, fs::Array{<:Int,2}, vals::::Array{<:Real,1}, level::Real)

Computes the "level curve" for a function defined as a collection of values on the vertices of a triangular mesh.
# Example

```julia-repl
julia> levelcurve(points,faces,values,level)
2×4 Array{Float64,2}:
 0.0  0.5  0.5  0.0
 0.5  1.0  1.0  0.5

julia> points = [[0 0];[1 0];[0 1];[1 1]]; faces = [[1 2 3];[2 4 3]]; values = [0,1,1,0]; level =1/2;

julia> levelcurve(points,faces,values,level)
2×4 Array{Float64,2}:
 0.0  0.5  0.5  0.0
 1.0  0.5  0.5  1.0
```

In the example above, each row of the matrix correspond to a level edge.
Columns 1 and 2  are the coordinates of the first vertex, and
columns 3 and 4 are the coordinates of the second vertex.
The edges are oriented so that the gradient of the function points to their left side.
"""
function levelcurve(ps::Array{<:Real,2}, fs::Array{<:Int,2}, vals::Array{<:Real,1}, level::Real)
    @assert size(fs,2) == 3
    dim = size(ps,2)
    law =[[[3 1],[3 2]], [[2 3],[2 1]], [[1 3],[1 2]],[[1 2],[1 3]],[[2 1],[2 3]],[[3 2],[3 1]]]
    cnt = 0
    es = zeros(size(fs,1),2*dim) # make enough space to hold edges (edges as pairs of points in 3-space)
    for i in 1:size(fs,1)
        fvs = fs[i,:] # Face vertices
        vals_fvs= vals[fvs] .> level
        ss = 1*vals_fvs[3] + 2*vals_fvs[2] + 4*vals_fvs[1]
        if 0 < ss < 7 # When intersection occurs
            ea, eb = law[ss]
            a1, a2 = fvs[ea]
            b1, b2 = fvs[eb]
            val_a1 = vals[a1]; val_a2 = vals[a2]
            val_b1 = vals[b1]; val_b2 = vals[b2]
            λa = (level-val_a1)/(val_a2-val_a1)
            λb = (level-val_b1)/(val_b2-val_b1)
            pa = (1-λa) *ps[a1,:] + λa*ps[a2,:]
            pb = (1-λb) *ps[b1,:] + λb*ps[b2,:]
            es[cnt+1,1:dim] = pa
            es[cnt+1,dim+1:2*dim] = pb
            cnt = cnt + 1
        end
    end
    es[1:cnt,:]
end

"""
    levelsurface(ps::Array{<:Real,2}, ts::Array{<:Int,2}, vals::Array{<:Real,1}, level::Real)

Computes the "level surface" for a function defined as a collection of values on the vertices of a tetrahedral mesh.
# Example

```julia-repl
julia> points = [[0 0 0];[1 0 0];[0 1 0];[0 0 1]]; tetrahedrons= [[1 2 3 4];]; values = [0,1,1,1]; level =1/2;

julia> levelsurface(points, tetrahedrons, val
valtype values
julia> levelsurface(points, tetrahedrons, values, level)
1×9 Array{Float64,2}:
 0.5  0.0  0.0  0.0  0.5  0.0  0.0  0.0  0.5
```

In the example above, each row of the matrix correspond to a level triangle.
Columns 1, 2 and 3 are the coordinates of the first vertex,
columns 4, 5 and 6 are the coordinates of the second vertex, and
columns 7, 8 and 9 are the coordinates of the third vertex.
The triangles are oriented so that the gradient of the function points to their left side.
"""
function levelsurface(ps::Array{<:Real,2}, ts::Array{<:Int,2}, vals::Array{<:Real,1}, level::Real)
    @assert size(ts,2) == 4
    dim = size(ps,2)
    # vertex1+
    l08 = [[[1 2],[1 4],[1 3]]]
    l07 = [[[1 2],[1 3],[1 4]]]
    # vertex2+
    l04 = [[[2 1],[2 3],[2 4]]]
    l11 = [[[2 1],[2 4],[2 3]]]
    # vertex3+
    l02 = [[[3 1],[3 4],[3 2]]]
    l13 = [[[3 1],[3 2],[3 4]]]
    # vertex4+
    l01 = [[[4 1],[4 2],[4 3]]]
    l14 = [[[4 1],[4 3],[4 2]]]

    #vertices??+
    l12 = [[[1 3],[2 3],[2 4]], [[1 3],[2 4],[1 4]]]
    l03 = [[[1 3],[2 4],[2 3]], [[1 3],[1 4],[2 4]]]
    #vertices??+
    l10 = [[[1 2],[1 4],[4 3]], [[1 2],[4 3],[2 3]]]
    l05 = [[[1 2],[4 3],[1 4]], [[1 2],[2 3],[4 3]]]
    #vertices??+
    l09 = [[[1 2],[2 4],[1 3]], [[1 3],[2 4],[4 3]]]
    l06 = [[[1 2],[1 3],[2 4]], [[1 3],[4 3],[2 4]]]

    ll = [l01,l02,l03,l04,l05,l06,l07,l08,l09,l10,l11,l12,l13,l14]
    cnt = 0
    fs = zeros(size(ts,1),3*dim) # make enough space to hold edges (edges as triplets of points in 3-space)

    for i in 1:size(ts,1)
        tvs = ts[i,:] # Tetrahedron vertices
        vals_tvs= vals[tvs] .> level
        c0, c1, c2, c3 = vals[tvs] .> level
        ss = 1*vals_tvs[4] + 2*vals_tvs[3] + 4*vals_tvs[2] + 8*vals_tvs[1]
        if 0 < ss < 15 # When intersection occurs
            for tr in ll[ss]
                ea, eb, ec = tr # edges having intersection points forming a level face
                a1, a2 = tvs[ea]
                b1, b2 = tvs[eb]
                c1, c2 = tvs[ec]
                val_a1 = vals[a1]; val_a2 = vals[a2]
                val_b1 = vals[b1]; val_b2 = vals[b2]
                val_c1 = vals[c1]; val_c2 = vals[c2]
                λa = (level-val_a1)/(val_a2-val_a1)
                λb = (level-val_b1)/(val_b2-val_b1)
                λc = (level-val_c1)/(val_c2-val_c1)
                pa = (1-λa) *ps[a1,:] + λa*ps[a2,:]
                pb = (1-λb) *ps[b1,:] + λb*ps[b2,:]
                pc = (1-λc) *ps[c1,:] + λc*ps[c2,:]
                fs[cnt+1,1:dim] = pa
                fs[cnt+1,dim+1:2*dim] = pb
                fs[cnt+1,2*dim+1:3*dim] = pc
                cnt = cnt + 1
            end
        end
    end
    fs[1:cnt,:]
end
