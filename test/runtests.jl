using HWfuncapp
using Test

@testset "HWfuncapp.jl" begin
	# test that q1 with 15 nodes has an error smaller than 1e-9
    	@test q1(15)[:error] < 1e-9
    # test that the integral of function h [-10,0] ≈ 1.46039789878568

    @test false  # by default fail
end
