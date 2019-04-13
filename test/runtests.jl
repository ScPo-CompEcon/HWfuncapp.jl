using HWfuncapp
using Test

@testset "HWfuncapp.jl" begin
	# test that q1 with 15 nodes has an error smaller than 1e-9
		@testset "test error" begin
			err = q1()
			@test err[:error] < 1e-9
		end
    # test that the integral of function h [-10,0] ≈ 1.46039789878568
		@testset "test integral" begin
			integral = q3(10)[2]
			@test integral ≈ 1.46039789878568
		end

    #@test false  # by default fail
end
