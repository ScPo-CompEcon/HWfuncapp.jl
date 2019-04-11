
module HWfuncapp

	using FastGaussQuadrature # to get chebyshevnodes

	# you dont' have to use PyPlot. I did much of it in PyPlot, hence
	# you will see me often qualify code with `PyPlot.plot` or so.
	using PyPlot
	using LinearAlgebra
	import ApproXD: getBasis, BSpline
	using Distributions
	using ApproxFun
	using Plots

	export scalemap, unitmap, TT, ChebyT, q1, q2, q3, q7, runall

	ChebyT(x,deg) = cos(acos(x)*deg)
	unitmap(x,lb,ub) = 2 .* (x .- lb) ./ (ub .- lb) .- 1	#[a,b] -> [-1,1]
	scalemap(x,lb,ub) = 0.5 .* (lb .+ ub) .+ 0.5 .* (ub .- lb) .* x #scaling to [a,b]

	function TT(x,J,lb,ub)
		M = length(x)
		Phi = zeros(M,J)
		for i in 1:M
			for j in 1:J
				Phi[i,j] = ChebyT(unitmap(x[i],lb,ub),j-1)
			end
		end
		return Phi
	end


	#function IPmatCh(n) #Chebyshev interpolation matrix
	#	M = zeros(n,n)
	#	for i in 1:n
	#		for j in 1:n
	#			M[i,j] = cos(((n - i + 0.5)*(j - 1)*(pi))/(n))
	#		end
	#	end
	#end

	#function cHat(n)
	#	return inv(IPmatCh(n)) * f(gch(-3,3,n))
	#end
	#f(x) = x .+ 2x.^2 - exp.(-x)

	#function gch(lb, ub, nnodes)
	#	return scalemap(gausschebyshev(nnodes)[1], lb, ub)
	#end


	function q1(n=15)
		lb = -3
		ub = 3
		n_new = 100
		f(x) = x .+ 2x.^2 - exp.(-x)
		chNodes = scalemap(gausschebyshev(n)[1], lb, ub)
		cIP = inv(TT(chNodes,n,lb,ub))*f(chNodes)
		grid = Vector(range(lb,length = n_new,stop = ub))
		true_new = f(grid)
		IP_new = (TT(grid,n,lb,ub))*cIP
		err = true_new .- IP_new
		p = Plots.plot(layout = 2, dpi = 400)
		Plots.plot!(p[1],grid,true_new,label = "True",lw = 1,linecolor = "black")
		Plots.plot!(p[1],grid,IP_new, label = "IP Approx", lw = 2,linestyle = :dot, linecolor = "red")
		Plots.plot!(p[2],grid,err, label = "error", lw = 1, linecolor = "green")
		# without using PyPlot, just erase the `PyPlot.` part
		Plots.savefig(p,joinpath(dirname(@__FILE__),"..","q1.png"))
		return Dict(:error=>maximum(abs,err))
	end

	function q2(b::Number = 4)
		@assert b > 0
		# use ApproxFun.jl to do the same:
		S = Chebyshev(-b..b)
		ft = Fun(x -> x .+ 2x.^2 - exp.(-x),-b..b)
		n,m = 100,50
		grid = Vector(range(-b, stop = b, length = n))
		true_v = ft.(grid)
		V = zeros(n,m)
		for k = 1:m
    		V[:,k] = Fun(S,[zeros(k-1);1]).(grid)
		end
		fa = Fun(S,V\true_v)
		approx_v = fa.(grid)
		err = approx_v .- true_v
		p = Plots.plot(layout = 2, dpi = 400)
		Plots.plot!(p[1],grid,true_v,label = "True",lw = 1,linecolor = "black")
		Plots.plot!(p[1],grid,approx_v, label = "IP Approx", lw = 2,linestyle = :dot, linecolor = "red")
		Plots.plot!(p[2],grid,err, label = "error", lw = 1, linecolor = "green")
		Plots.savefig(p,joinpath(dirname(@__FILE__),"..","q2.png"))
	end

	function q3(b::Number = 10)
		@assert b > 0
		x = Fun(identity,-b..b)
		f = sin(x^2)
		g = cos(x)
		h = f - g
		rh = roots(h)
		p = Plots.plot(h, label = "h(x)", dpi = 400, lw = 1, linecolor = "black")
		Plots.scatter!(p,rh,h.(rh), label = "roots of h(x)", markercolor = "green", markersize = 5)
		Plots.savefig(p,joinpath(dirname(@__FILE__),"..","q3.png"))
		y = Fun(identity,-10..0)
		hh = sin(y^2) - cos(y)
		hh_int1 = cumsum(hh)
		hh_int = hh_int1 + hh(-10)
		integral = ApproxFun.norm(hh_int(0) - hh_int(-10))
		return (p,integral)
	end

	# optinal
	function q4()

		return fig
	end


	# I found having those useful for q5
	mutable struct ChebyType
		f::Function # fuction to approximate
		nodes::Union{Vector,LinRange} # evaluation points
		basis::Matrix # basis evaluated at nodes
		coefs::Vector # estimated coefficients

		deg::Int 	# degree of chebypolynomial
		lb::Float64 # bounds
		ub::Float64

		# constructor
		function ChebyType(_nodes::Union{Vector,LinRange},_deg,_lb,_ub,_f::Function)
			n = length(_nodes)
			y = _f(_nodes)
			_basis = Float64[ChebyT(unitmap(_nodes[i],_lb,_ub),j) for i=1:n,j=0:_deg]
			_coefs = _basis \ y  # type `?\` to find out more about the backslash operator. depending the args given, it performs a different operation
			# create a ChebyComparer with those values
			new(_f,_nodes,_basis,_coefs,_deg,_lb,_ub)
		end
	end

	# function to predict points using info stored in ChebyType
	function predict(Ch::ChebyType,x_new)

		true_new = Ch.f(x_new)
		basis_new = Float64[ChebyT(unitmap(x_new[i],Ch.lb,Ch.ub),j) for i=1:length(x_new),j=0:Ch.deg]
		basis_nodes = Float64[ChebyT(unitmap(Ch.nodes[i],Ch.lb,Ch.ub),j) for i=1:length(Ch.nodes),j=0:Ch.deg]
		preds = basis_new * Ch.coefs
		preds_nodes = basis_nodes * Ch.coefs

		return Dict("x"=> x_new,"truth"=>true_new, "preds"=>preds, "preds_nodes" => preds_nodes)
	end

	function q5(deg=(5,9,15),lb=-1.0,ub=1.0)

		runge(x) = 1.0 ./ (1 .+ 25 .* x.^2)


		PyPlot.savefig(joinpath(dirname(@__FILE__),"..","q5.png"))

	end



	function q6()

		# compare 2 knot vectors with runge's function

		PyPlot.savefig(joinpath(dirname(@__FILE__),"..","q6.png"))

	end

	function q7(pr = 0.025)
		deg = 3
		lb = -1
		ub = 1
		nk = 13
		nev = 65
		f(x) = (abs.(x)).^0.5
		bs1 = BSpline(nk,deg,lb,ub)
		knots2 = vcat(range(lb,stop = -pr, length = 5),0,0,0,range(pr,stop = ub, length = 5))
		bs2 = BSpline(knots2,deg)
		eval = Vector(range(lb,stop = ub, length = nev))
		c1 = getBasis(eval,bs1) \ f(eval)
		c2 = getBasis(eval,bs2) \ f(eval)
		grid = Vector(range(lb, stop = ub, length = 1000))
		tr = f(grid)
		a1 = getBasis(grid, bs1) * c1
		a2 = getBasis(grid, bs2) * c2
		err1 = a1 .- tr
		err2 = a2 .- tr
		p = Plots.plot(layout = 3, dpi = 400)
		Plots.plot!(p[1],grid,tr,label = "True",lw = 1.5,linecolor = "black")
		Plots.plot!(p[2],grid,a1, label = "Approx1", lw = 1.5, linecolor = "blue")
		Plots.plot!(p[2],grid,a2, label = "Approx2", lw = 1.5, linecolor = "red")
		Plots.plot!(p[3],grid,err1, label = "error1", lw = 1.5, linecolor = "blue")
		Plots.plot!(p[3],grid,err2, label = "error2", lw = 1.5, linecolor = "red")
		Plots.savefig(p,joinpath(dirname(@__FILE__),"..","q7.png"))
		return (p,maximum(abs,err2))
	end

	#function q7bis()
	#	G = Vector(range(0.001, stop = 0.2, length = 1000))
	#	return Plots.plot(G,q7.(G))
	#end


	# function to run all questions
	function runall()
		@info("running all questions of HW-funcapprox:")
		q1(15)
		q2(3)
		q3(10)
		#q4()
		#q5()
		#q6()
		q7()
	end



end # module
