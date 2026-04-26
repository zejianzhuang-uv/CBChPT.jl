
using DataFrames
using DataStructures
using Statistics
include("./mB_nlo.jl")
include("./mB_one_loop.jl")
include("./error_analysis.jl")




function mB_up_to_p3_errband(par_sample::AbstractDataFrame, mphi::AbstractDataFrame)
    mB_sample = []
    for par in eachrow(par_sample)
        m0, b0, bD, bF = Tuple(par)
        LECs_dict = OrderedDict(:b0 => b0, :bD => bD, :bF => bF)
        mB = mB_up_to_p3(m0, LECs_dict, mphi)
        push!(mB_sample, mB)
    end
    mB_sample = stack(mB_sample)
    mB_err = std(mB_sample, dims=3)[:, :, 1] |> loess_smooth
    name = (size(mB_err, 2) == 5 ? [:mN_err, :mL_err, :mS_err, :mXi_err, :mBsu3_err] : [:mN_err, :mL_err, :mS_err, :mXi_err])
    df = DataFrame(mB_err, name)
    return df
end





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
    mB = zeros(Float64, n, ncol)

    mpi  = mphi.mpi
    mK   = mphi.mK
    meta = mphi.meta
    fpi  = mphi.fpi
    fK   = mphi.fK
    feta = mphi.feta
    fsu3 = has_su3 ? mphi.fsu3 : nothing

    for i in 1:n
        mB[i, 1:4] = mB_up_to_p3(m0, mpi[i], mK[i], meta[i], fpi[i], fK[i], feta[i], LECs_dict)
        has_su3 && (mB[i, 5] = mB_up_to_p3_su3(m0, mpi[i], fsu3[i], LECs_dict))
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
    b0, bD, bF = (LECs_dict[:b0]*1e-3, LECs_dict[:bD]*1e-3, LECs_dict[:bF]*1e-3) # change to MeV^-1
    mB2 = CBnlo_su3(b0, bD) * mpi^2
    mB2ol = CB_one_loop_su3 * H3(mpi, fpi, m0)
    mB = m0 + mB2 + mB2ol
    return mB
end