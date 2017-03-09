

module IntTest

	using HW_int
	using Base.Test

	@test 1==1

	using Base.Test

println("Running the Tests")

## TESTS on EXERCICE 1
######################
n = 10
p1 = 4
p2  = 1

## Testing the Demand Function

@testset "Demand function" begin
		 @test HW_int.demand(1) == 2
		 @test HW_int.demand(4) == 1
end

# Testing the Weights & Integral Approximation Functions

@testset "Weights Test" begin

	x = HW_int.question_1b(n)
	# For the Gauss-Legendre quadrature, the weights are defined for nodes in [-1;1] so we test that the weighted sum of a constant indeed returns the integral of this constant over [-1; 1] (i.e 2*constant)
	constant = 3
	@test dot(x[:2],ones(n)*constant) ≈ constant * 2

	# For the function question_1c, and question_1d we simply test that the sum of the weights equal to 1
	@test sum(HW_int.question_1c(n)[:2]) ≈ 1
	@test sum(HW_int.question_1d(n)[:2]) ≈ 1

end

## TESTS on Computing Time
##########################
println("Computing effectiveness:")
@time HW_int.question_1b(n)
@time HW_int.question_1c(n)
@time HW_int.question_1d(n)

println("End of tests")
println("--------------------")

end
