
module HW_int

	using FastGaussQuadrature
	using Roots
	using Sobol
	using Gadfly
	using Distributions

	# Question 1 parameters
	p1 = 4
	p2 = 1

	# demand function
	function demand(p)
	  return 2*(p)^(-1/2)
	end

	# Question 2 paramters
	Supply = 2
	m = [0; 0]
	varcov = [0.02 0.01; 0.01 0.01]


################ QUESTION 1 #######################

	function question_1b(n)

		nodes0 = gausslegendre(n)[:1]
		weights = gausslegendre(n)[:2]
		nodes = nodes0*(p1-p2)/2 + (p1+p2)/2 # Adjusted nodes : from the [-1;1] space to the [1;4] space

		# Transformation of the demand function
		result = (p1 - p2)/2 * dot(weights, demand.(nodes))
		println("Gauss Legendre Quadrature solution : ", result)

		global plotb = Gadfly.plot(layer(x=nodes, y=demand.(nodes), Geom.point, Theme(default_point_size = 2pt, default_color=colorant"blue")),
				layer(demand,0,4), Guide.xlabel("Prices"), Guide.ylabel("Demand"),Guide.Title("Question 1b : Gauss-Legendre Quadrature"))
		######## We displayed the figures outside the question_1b function because
		######## otherwise it would print 2 or 3 times the same graph

		return nodes, weights

	end


	function question_1c(n)

		# Drawing nodes randomly from a uniform law
		#srand(123)
		nodes = rand(n)*(p1-p2)+p2
	  weights = collect(linspace(1/n,1/n,n))

	  result = (p1-p2)*dot(weights,demand.(nodes)) #Nb : the result is multiplied by (p1-p2) because of the integration by substitution we do : from a variable in [0;1] to a variable in [p2;p1]

		println("Monte Carlo solution : ", result)
		global plotc = Gadfly.plot(layer(x=nodes, y=demand.(nodes), Geom.point, Theme(default_point_size = 2pt, default_color=colorant"blue")),
		  layer(demand,0,4), Guide.xlabel("Prices"), Guide.ylabel("Demand"),Guide.Title("Question 1c : Monte Carlo Simulation"))
		######## We displayed the figures outside the question_1b function because
		######## otherwise it would print 2 or 3 times the same graph

		return nodes, weights

	end


	function question_1d(n)

		# Generating the Sobol sequence grid
	  nodes=zeros(n)
	  s =  SobolSeq(1,p2,p1)
	  for i in 1:n
	    nodes[i] = next(s)[1]
	  end
	  weights = collect(linspace(1/n,1/n,n))

	  result = (p1-p2)*dot(weights, demand.(nodes)) #Nb : the result is multiplied by (p1-p2) because of the integration by substitution we do : from a variable in [0;1] to a variable in [p2;p1]
	  println("Pseudo Monte-Carlo solution : ", result)

		global plotd = Gadfly.plot(layer(x=nodes, y=demand.(nodes), Geom.point, Theme(default_point_size=2pt,default_color=colorant"blue")),
		  layer(demand,0,4),  Guide.xlabel("Prices"), Guide.ylabel("Demand"), Guide.Title("Question 1d : Pseudo Monte Carlo Simulation"))
		######## We displayed the figures outside the question_1b function because
		######## otherwise it would print 2 or 3 times the same graph

		return nodes, weights

	end


################ QUESTION 2 #####################

	function question_2a(n)

		# 1) Approximation of the Expect Demand functions (domestic and exportation) via Gauss-Hermite quadrature

		# Defining a 2 dimension set of nodes and associated weights
		nodes0 = Matrix(1,2)
		nodes0 = hcat(repeat(gausshermite(n)[:1],inner=[1],outer=[n]), repeat(gausshermite(n)[:1],inner=[n],outer=[1]))
		weights = kron(gausshermite(n)[:2],gausshermite(n)[:2])./pi

		# Adapting the nodes to take into account the correlation between the 2 shocks
		omega = cholfact(varcov)[:L]
		grid = nodes0 *omega + repmat(transpose(m),n^2)

		# Gauss-Hermite Expectation of expected total surplus
		# (Actually what we simulated is log(theta) which follows a normal distribution - and not theta.So we express the demand in terms of log(theta) and not theta)
		Expected_Surplus = (p) -> dot(weights,exp(grid[:,1] - log(p))) +  dot(weights,exp(grid[:,2] - log(p))) - Supply

		# 2) Equilibrium Price Determination (Expected demand = Supply) and Variance
		Expected_p = fzero(x->Expected_Surplus(x),0.05,5)
		Var_p= fzero(x -> Expected_Surplus(x^2)-Expected_p^2,0.05,5)

		println("Expected Equilibrium Price (Gauss-Hermite method): ", Expected_p)
		println("Variance of the Equilibium Price:", Var_p)
		return grid

	end


	function question_2b(n)


		# 1) Monte Carlo Simulation of the expected surplus (X+D)

		# Definning a 2-dimension grid of N^2 points
		srand(323)
		nodes0 = hcat(kron(ones(n),randn(n)),kron(randn(n),ones(n)))
		weights = kron(collect(linspace(1/n,1/n,n)),collect(linspace(1/n,1/n,n)))

		# Adapting the nodes to take into account the correlation between the 2 shocks
		omega = cholfact(varcov)[:L]
		grid = nodes0 *omega + repmat(transpose(m),n^2)

		# Expected surplus
	  Expected_Surplus = (p) -> dot(weights,exp(grid[:,1] -log(p))) + dot(weights,exp(grid[:,1] -log(p))) - Supply

		# 2) Equilibrium Price Determination (Expected demand = Supply)
		Expected_p = fzero(Expected_Surplus,0.05,5)
		Var_p= fzero(x -> Expected_Surplus(x^2)-Expected_p^2,0.05,5)

		println("Expected Equilibrium Price (Monte-Carlo method) : ", Expected_p)
		println("Variance of the Equilibium Price (Monte-Carlo method):", Var_p)
#		plot()

	end


	# function to run all questions
	function runall(n=10)

		println("running all questions of HW-integration:")

		println("")
		println("Results of question 1:")
		println("----------------------")
		question_1b(n)
		display(plotb)
		question_1c(n)
		display(plotc)
		question_1d(n)
		display(plotd)

		println("")
		println("Results of question 2:")
		println("----------------------")
		question_2a(n)
		question_2b(n)
		#q2 = question_2a(n)
		#println(q2)
		#q2b = question_2b(n)
		#println(q2b)
		#println("end of HW-integration")

	end

end
