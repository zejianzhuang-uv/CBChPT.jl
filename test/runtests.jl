using CBChPT
using DataStructures
using Test

@testset "CBChPT.jl" begin
    m0, b0, bD, bF = 805.0409408126626, -0.664609603537589, 0.062434962905761365, -0.35384103271852047
    LECs_dict = OrderedDict(
        :b0 => b0,
        :bD => bD,
        :bF => bF
    )
    mpi, mK, meta = 138., 495., 547.
    fpi, fK, feta = 93., 108., 120.
    
    println((LECs_dict))
    mB_up_to_p3(m0, mpi, mK, meta, fpi, fK, feta, LECs_dict) |> println
    println(mB_up_to_p3_su3(m0, mpi, fpi, LECs_dict))
end
