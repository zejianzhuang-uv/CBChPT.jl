

"""
(p^3)
"""
function mB_one_loop(mphi::Vector{Float64}, fphi::Vector{Float64}, mu::Float64)
    h = @. H3(mphi, fphi, mu)
    return CB_one_loop * h
end





function H3(mϕ::Float64, fϕ::Float64, mu::Float64)
    t1 = -2 * mϕ^3 * (sqrt(1 - mϕ^2 / (4mu^2) ) * acos(mϕ / (2mu) ) + mϕ/(2mu) * log(mϕ / mu) )
    return 1 / (4π * fϕ)^2 * t1
end




const D = 0.8
const F = 0.46

CB_one_loop = zeros(Float64, (4, 3))
CB_one_loop[1, :] = [3/2(D+F)^2, 1/3*(5D^2-6*D*F+9F^2), 1/6*(D-3F)^2] # N
CB_one_loop[2, :] = [2D^2, 2/3*(D^2+9F^2), 2/3*D^2]
CB_one_loop[3, :] = [2/3*(D^2+6F^2), 2*(D^2+F^2), 2/3*D^2]
CB_one_loop[4, :] = [3/2*(D-F)^2, 1/3*(5D^2+6*D*F+9*F^2), 1/6*(D+3F)^2]

