# Utility functions here

"""
    get_coords(points::Array{<:Real,2})

Gets the coordinates of a colleciton of points stored in the rows of a matrix.

# Example
```julia-repl
julia> points = [[1.0 2.0];[9.4 1.7]];

julia> x,y = get_coords(points);

julia> print(x)
[1.0, 9.4]
julia> print(y)
[2.0, 1.7]
```
"""
function get_coords(points::Array{<:Real,2})
    [points[:,i] for i in 1:size(points,2)]
end

function get_points(coords)
    [coords[i,:] for i=1:size(coords,1)]
end

function get_logistic(;x0=0,L=1,k=1)
    function out(x)
        return L/(1+exp(-k*(x-x0)))
    end
    return out
end
