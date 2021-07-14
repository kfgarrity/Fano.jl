using Fano
using Test

function do_test()

    x = Fano.main_fit_VV()
    @test x[1] ≈ 0.6742588756733324 * sqrt(200)
end

do_test()
