######
#
# 2-dimensional stuff.
#

function regulargrid2d(box, res)
    xmin, ymin, xmax, ymax = box
    rx, ry = res
    dx = (xmax-xmin)/(rx-1)
    dy = (ymax-ymin)/(ry-1)

    vs_cnt = rx*ry # vertices count
    es_cnt = (rx-1)*ry + rx*(ry-1)+ (rx-1)*(rx-1) # Horizontal + vertical + diagonal
    fs_cnt = 2*(rx-1)*(ry-1) # down + up

    ps = zeros(vs_cnt, 3) # points
    vs = zeros(Int,vs_cnt,1) # vertices
    es = zeros(Int,es_cnt,2) # edges
    fs = zeros(Int,fs_cnt,3) # faces

    # rows, columns and values for edge-vertex adjacency
    e_rows = zeros(Int, 2*es_cnt)
    e_cols = zeros(Int, 2*es_cnt)
    e_vals = zeros(Int, 2*es_cnt)

    # rows, columns and values for face-edge adjacency
    f_rows = zeros(Int, 3*fs_cnt)
    f_cols = zeros(Int, 3*fs_cnt)
    f_vals = zeros(Int, 3*fs_cnt)

    esc = 1
    fsc = 1
    for r in 0:(ry-1)
        for c in 0:(rx-1)
            # vertices indices
            v00 = r*rx+ c
            v01 = v00 + 1
            v10 = v00 + rx
            v11 = v01 + rx
            # horizontal edges
            he00 = r*(rx-1) + c
            he01 = he00 + 1
            he10 = he00 + (rx-1)
            he11 = he01 + (rx-1)
            # vertical edges
            ve00 = (rx-1)*ry + r*rx + c
            ve01 = ve00 + 1
            ve10 = ve00 + rx
            ve11 = ve01 + rx
            # diagonal edges
            de00 = (rx-1)*ry + rx*(ry-1) + r*(rx-1) + c
            de01 = de00 + 1
            de10 = de00 + (rx-1)
            de11 = de01 + (rx-1)
            # down faces
            df00 = r*(rx-1) + c
            df01 = df00 + 1
            df10 = df00 + (rx-1)
            df11 = df01 + (rx-1)
            # up faces
            uf00 = (rx-1)*(ry-1) + r*(rx-1) + c
            uf01 = uf00 + 1
            uf10 = uf00 + (rx-1)
            uf11 = uf01 + (rx-1)

            # setting points
            ps[v00+1,:] = [xmin + c*dx, ymin + r*dy, 0.0]
            # setting vertices
            vs[v00+1] = v00
            # Setting edges and faces
            if c < (rx-1) # horizontal edges
                es[he00+1,:] = [v00,v01]
                e_rows[esc] = he00; e_cols[esc] = v00; e_vals[esc] = -1; esc = esc + 1
                e_rows[esc] = he00; e_cols[esc] = v01; e_vals[esc] = 1; esc = esc + 1
            end
            if r < (ry-1) # vertical edges
                es[ve00+1,:] = [v00,v10]
                e_rows[esc] = ve00; e_cols[esc] = v00; e_vals[esc] = -1; esc = esc + 1
                e_rows[esc] = ve00; e_cols[esc] = v10; e_vals[esc] = 1; esc = esc + 1
            end
            if r <(ry-1) && c < (rx-1)
                fs[df00+1,:] = [v00,v01,v11] # down faces
                f_rows[fsc] = df00; f_cols[fsc] = he00; f_vals[fsc] = 1; fsc = fsc + 1
                f_rows[fsc] = df00; f_cols[fsc] = ve01; f_vals[fsc] = 1; fsc = fsc + 1
                f_rows[fsc] = df00; f_cols[fsc] = de00; f_vals[fsc] = -1; fsc = fsc + 1

                fs[uf00+1,:] = [v11,v10,v00] # up faces
                f_rows[fsc] = uf00; f_cols[fsc] = de00; f_vals[fsc] = 1; fsc = fsc + 1
                f_rows[fsc] = uf00; f_cols[fsc] = he10; f_vals[fsc] = -1; fsc = fsc + 1
                f_rows[fsc] = uf00; f_cols[fsc] = ve00; f_vals[fsc] = -1; fsc = fsc + 1

                es[de00+1,:] = [v00, v11]    # diagonal edges
                e_rows[esc] = de00; e_cols[esc] = v00; e_vals[esc] = -1; esc = esc + 1
                e_rows[esc] = de00; e_cols[esc] = v11; e_vals[esc] = 1; esc = esc + 1
            end
        end
    end
    ps, vs.+1, es.+1, fs.+1
end

######
#
# 3-dimensional stuff.
#
function regulargrid3d(box, res)
    xmin, ymin, zmin, xmax, ymax, zmax = box
    resx, resy, resz = res
    resxy = resx*resy
    resxyz = resxy*resz
    dx = (xmax-xmin)/(resx-1)
    dy = (ymax-ymin)/(resy-1)
    dz = (zmax-zmin)/(resz-1)
    ps = zeros(resxyz, 3)
    vs = zeros(Int,resxyz,1)
    ts = zeros(Int,6*(resx-1)*(resy-1)*(resz-1),4)
    tc = 1
    for k in 0:(resz-1)
        for j in 0:(resy-1)
            for i in 0:(resx-1)
                id = resxy*k+ resx*j + i
                vs[id+1] = id
                ps[id+1,:] = [xmin+i*dx,ymin+j*dy,zmin+k*dz]
                v0 = resxy*k + j*resx + i
                v1 = resxy*k + j*resx + i + 1
                v2 = resxy*k + (j+1)*resx + i+1
                v3 = resxy*k + (j+1)*resx + i

                v4 = resxy*(k+1) + j*resx + i
                v5 = resxy*(k+1) + j*resx + i + 1
                v6 = resxy*(k+1) + (j+1)*resx + i + 1
                v7 = resxy*(k+1) + (j+1)*resx + i
                if i <(resx-1) && j < (resy-1) && k < (resz-1)
                    ts[tc+0,:] = [v2,v0,v3,v7]
                    ts[tc+1,:] = [v0,v2,v6,v7]
                    ts[tc+2,:] = [v4,v0,v6,v7]
                    ts[tc+3,:] = [v6,v0,v1,v2]
                    ts[tc+4,:] = [v0,v6,v1,v4]
                    ts[tc+5,:] = [v6,v5,v1,v4]
                    tc = tc + 6
                end
            end
        end
    end
    ps, ts.+1
end

function compute_subfaces(faces::Array{Int,2})
    fs_cnt, k = size(faces)
    mask = ones(Bool,k)
    sfs = Array{Array{Int64,2}}(undef,k*fs_cnt)
    for i in 1:fs_cnt
        for j in 1:k
            mask[j] = false
            sfs_index = (i-1)*k + j
            sfs[sfs_index] = sort(faces[[i],mask],dims=2)
            mask[j] = true
        end
    end
    subfaces = vcat(unique(sfs)...)
    subfaces
end
