
using Loess



function loess_smooth(f::AbstractVector{Float64})
    n = length(f)
    predict(loess(Float64.(1:n), f), Float64.(1:n))
end

function loess_smooth(f::AbstractMatrix{Float64})
    mapslices(loess_smooth, f, dims=1)
end