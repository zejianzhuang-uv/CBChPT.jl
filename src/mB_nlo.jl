




function mBnlo(mpiKsq::AbstractVector{Float64}, b0::Float64, bD::Float64, bF::Float64; su3=false)
    if su3 == false
        C = CBnlo(b0, bD, bF)
        # ein"ij,i->j"(C, mpiKsq)
        C*mpiKsq
    end
end


"""
b0, bD, bF have the unit GeV^-1.
b0 = -0.665, bD = 0.062, bF = -0.354
"""
function CBnlo(b0::Float64, bD::Float64, bF::Float64)
    C = zeros(Float64, (4, 2))
    C[1, :] = - [2b0 + 4bF, 4b0+4bD-4bF]
    C[2, :] = - [2/3*(3b0-2bD), 2/3*(6b0+8bD)]
    C[3, :] = - [(2b0+4bD), 4b0]
    C[4, :] = - [(2b0-4bF), (4b0+4bD+4bF)]
    return C
end

function CBnlo_su3(b0::Float64, bD::Float64)
    return -2(3b0+2bD)
end