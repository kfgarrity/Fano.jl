using Fano
using Test

function do_test()

    x = Fano.main_fit_VV()
    @test x[1] ≈ 40.56858406766799
end

do_test()
