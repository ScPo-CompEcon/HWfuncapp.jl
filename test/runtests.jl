using HWfuncapp
using Test

@testset "HWfuncapp.jl" begin
	# test that q1 with 15 nodes has an error smaller than 1e-9

    # test that the integral of function h [-10,0] ≈ 1.46039789878568
    @test q1(15)[:error]<1e-9
    @test (q3(10)[2] ≈ 1.46039789878568)== false # by default fail
	@test q2(4)[:error]<1e-9
end
