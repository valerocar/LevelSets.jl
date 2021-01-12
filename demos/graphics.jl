using GeometryBasics
using GeometryBasics: Mesh
using Colors
#
# Meshes in 3d space
#

"Compute a mesh from points (vertices) and faces"
function mesh3d(ps,fs; vertexColors = nothing )
    pps = [Point(ps[i,1],ps[i,2],ps[i,3]) for i in 1:size(ps)[1]]
    ffs = [TriangleFace(fs[i,1],fs[i,2],fs[i,3]) for i in 1:size(fs)[1]]
    if(vertexColors != nothing)
        return  meta(Mesh(pps, ffs),vertexColors=vertexColors)
    end
    Mesh(pps, ffs)
end

"Compute Mesh from the coordinates of triangles"
function mesh3d(faces_coords::Array{<:Real,2})
    nps = Point{3, Float64}[]
    nfs = NgonFace{3,Int}[]
    normals = Point{3,Float64}[]
    vcnt = 1
    for i in 1:size(faces_coords,1)
        p1 = faces_coords[i,1:3]
        p2 = faces_coords[i,4:6]
        p3 = faces_coords[i,7:9]
        push!(nps, Point(p1[1],p1[2],p1[3]))
        push!(nps, Point(p2[1],p2[2],p2[3]))
        push!(nps, Point(p3[1],p3[2],p3[3]))
        push!(nfs, TriangleFace(vcnt,vcnt+1,vcnt+2))
        vcnt = vcnt + 3
    end
    Mesh(nps,nfs)
end

"Compute line segments from vertices coordinates and edges"
function edges3d(ps,es; color=RGBA(1,1,1,1), linewidth=1)
    xs = ps[:,1]
    ys = ps[:,2]
    zs = ps[:,3]
    xes = transpose(xs[es])
    yes = transpose(ys[es])
    zes = transpose(zs[es])
    xes = vec(xes)
    yes = vec(yes)
    zes = vec(zes)
    MeshCat.LineSegments(Point.(xes,yes,zes),LineBasicMaterial(color=color, linewidth=linewidth))
end

"Compute line segments from the edges coordinates"
function edges3d(es_ps::Array{<:Real,2}; color=RGBA(1,1,1,1), linewidth=1)
    cnt = 2*size(es_ps,1)
    x = zeros(cnt)
    y = zeros(cnt)
    z = zeros(cnt)
    x[1:2:cnt] = es_ps[:,1]
    x[2:2:cnt] = es_ps[:,4]
    y[1:2:cnt] = es_ps[:,2]
    y[2:2:cnt] = es_ps[:,5]
    z[1:2:cnt] = es_ps[:,3]
    z[2:2:cnt] = es_ps[:,6]
    mat = LineBasicMaterial(color=color, linewidth=linewidth)
    LineSegments(Point.(x,y,z),mat)
end
