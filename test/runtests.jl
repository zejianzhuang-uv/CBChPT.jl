using CBChPT
using DataStructures
using Test

@testset "CBChPT.jl" begin
    LECs_dict = OrderedDict(
        :b0 => -0.665,
        :bD => 0.062,
        :bF => -0.354
    )
    mpi, mK, meta = 138., 495., 547.
    fpi, fK, feta = 93., 108., 120.
    m0 = 805e0
    println((LECs_dict))
    mB_up_to_p3(m0, mpi, mK, meta, fpi, fK, feta, LECs_dict) |> println
end
