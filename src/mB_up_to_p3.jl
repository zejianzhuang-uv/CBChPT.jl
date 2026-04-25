
using DataFrames

include("./mB_nlo.jl")
include("./mB_one_loop.jl")


function df_mB_up_to_p3(m0::Float64, LECs_dict::AbstractDict, mphi::AbstractDataFrame)
    mB = mB_up_to_p3(m0::Float64, LECs_dict::AbstractDict, mphi::AbstractDataFrame)
    name = (size(mB, 2) ? [:mN, :mL, :mS, :mXi] : [:mN, :mL, :mS, :mXi, :mBsu3])
    df = DataFrame(mB, name)
    insertcols!(df, 1, :mpi => mphi.mpi)
    return df
end

function mB_up_to_p3(m0::Float64, LECs_dict::AbstractDict, mphi::AbstractDataFrame)
    n = size(mphi, 1)
    has_su3 = hasproperty(mphi, :fsu3)
    ncol = has_su3 ? 5 : 4
    mB = zeros(Float64, (n, ncol))
    for i in 1:n
        mB[i, 1:4] = mB_up_to_p3(m0, mphi.mpi[i], mphi.mK[i], mphi.meta[i], mphi.fpi[i], mphi.fK[i], mphi.feta[i], LECs_dict)
        if has_su3
            mB[i, 5] = mB_up_to_p3_su3(m0, mphi.mpi[i], mphi.fsu3[i], LECs_dict)
        end
    end
    return mB
end


function mB_up_to_p3(m0::Float64, mpi::Float64, mK::Float64, meta::Float64, fpi::Float64, fK::Float64, feta::Float64, LECs_dict::AbstractDict)
    b0, bD, bF = (LECs_dict[:b0]*1e-3, LECs_dict[:bD]*1e-3, LECs_dict[:bF]*1e-3) 
    mB2 = mBnlo(Float64[mpi^2, mK^2], b0, bD, bF)
    mBol = mB_one_loop(Float64[mpi, mK, meta], Float64[fpi, fK, feta], m0)
    mB = @. m0 + mB2 + mBol
    return mB
end

function mB_up_to_p3_su3(m0::Float64, mpi::Float64, fpi::Float64, LECs_dict::AbstractDict)
    b0, bD, bF = (LECs_dict[:b0]*1e-3, LECs_dict[:bD]*1e-3, LECs_dict[:bF]*1e-3) 
    mB2 = CBnlo_su3(b0, bD) * mpi^2
    mB2ol = CB_one_loop_su3 * H3(mpi, fpi, m0)
    mB = m0 + mB2 + mB2ol
    return mB
end