module HWfuncapp

using FastGaussQuadrature # to get chebyshevnodes

# you dont' have to use PyPlot. I did much of it in PyPlot, hence
# you will see me often qualify code with `PyPlot.plot` or so.
using PyPlot  
import ApproXD: getBasis, BSpline
using LinearAlgebra
using Distributions
using ApproxFun
using Plots

ChebyT(x,deg) = cos(acos(x)*deg)
unitmap(x,lb,ub) = 2 .* (x .- lb) ./ (ub .- lb) .- 1	#[a,b] -> [-1,1]

export q1, q2, q3, q4, q5, q7, runall

function q1(n=15)
	f(x) = x .+ 2x.^2 - exp.(-x)

	deg = n:-1:1
	points = range(-3, 3, length=n)

	node = 3 * cos.((2 .* deg .- 1) .* pi ./ (2 .* n))
	y = f(node)

	cheb = zeros(n, n)
	for i in 1:n
		cheb[:, i] = ChebyT.(unitmap(node, -3, 3), i - 1)
	end
	c = cheb^-1 * y

	estim = zeros(n)
	for i in 1:n
		estim += c[i] * ChebyT.(unitmap(points, -3, 3), i - 1)
	end

	subplot(122)
	PyPlot.plot(points, f.(points))
	PyPlot.scatter(points, estim, color="red", s=2)
	xlabel("x")
	ylabel("y")
	subplot(121)
	PyPlot.plot(points, f.(points) - estim)
	# without using PyPlot, just erase the `PyPlot.` part
	PyPlot.savefig(joinpath(dirname(@__FILE__),"..","q1.png"))

	err = sum(abs.(f.(points) - estim))
	return Dict(:error=>maximum(abs,err))
end

function q2(b::Number)
	@assert b > 0
	f(x) = x .+ 2x.^2 - exp.(-x)

	precision = 100
	points = range(-3, 3, length=precision)
	# use ApproxFun.jl to do the same:
	x = Fun(f, Chebyshev(-b..b))
	estim = x.(points)

	subplot(122)
	PyPlot.plot(points, f.(points))
	PyPlot.scatter(points, estim, color="red", s=2)
	xlabel("x")
	ylabel("y")
	subplot(121)
	PyPlot.plot(points, f.(points) - estim)

	PyPlot.savefig(joinpath(dirname(@__FILE__),"..","q2.png"))
end

function q3(b::Number)
	x = Fun(identity, -b..b)
	f = sin(x^2)
	g = cos(x)
	h = f - g

	r = roots(h)

	p = Plots.plot(h, labels="h")
	Plots.scatter!(r, h.(r), labels="roots.h")
	Plots.savefig(joinpath(dirname(@__FILE__),"..","q3.png"))
	# p is your plot
	g = cumsum(h)
	integral = g(0) - g(-b)
	return (p,integral)
end

# optinal
function q4()
	precision = 10000
	points = range(-1, stop = 1, length=precision)
	fig = PyPlot.plot(points, ChebyT.(points, 0), label=0)
	for i in 1:8
		fig = PyPlot.plot(points, ChebyT.(points, i), label=i)
	end
	title("Chebyshev basis")
	legend()
	PyPlot.savefig(joinpath(dirname(@__FILE__),"..","q4.png"))
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
		_basis = Float64[ChebyT(unitmap(_nodes[i],_lb,_ub),j) for i=1:n,j=0:_deg - 1]
		_coefs = _basis \ y  # type `?\` to find out more about the backslash operator. depending the args given, it performs a different operation
		# create a ChebyComparer with those values
		new(_f,_nodes,_basis,_coefs,_deg,_lb,_ub)
	end
end

# function to predict points using info stored in ChebyType
function predict(Ch::ChebyType,x_new)

	true_new = Ch.f(x_new)
	basis_new = Float64[ChebyT(unitmap(x_new[i],Ch.lb,Ch.ub),j) for i=1:length(x_new),j=0:Ch.deg - 1]
	basis_nodes = Float64[ChebyT(unitmap(Ch.nodes[i],Ch.lb,Ch.ub),j) for i=1:length(Ch.nodes),j=0:Ch.deg - 1]
	preds = basis_new * Ch.coefs
	preds_nodes = basis_nodes * Ch.coefs

	return Dict("x"=> x_new,"truth"=>true_new, "preds"=>preds, "preds_nodes" => preds_nodes)
end

function q5(deg=(5,9,15),lb=-1.0,ub=1.0)

	runge(x) = 1.0 ./ (1 .+ 25 .* x.^2)
	f(x) = x .+ 2x.^2 - exp.(-x)
	precision = 10000
	points = range(lb, ub, length=precision)
	cheb_m = zeros(precision, length(deg))
	uni_m = zeros(precision, length(deg))

	for (i, d) in enumerate(deg)
		degse = d:-1:1
		node = 1/2 .* (ub .+ lb) .+ 1/2 .* (ub - lb) .* cos.((2 .* degse .- 1) .* pi ./ (2 .* d))
		uni = LinRange(lb, ub, d)
		cheb = ChebyType(node, d, lb, ub, runge)
		unic = ChebyType(uni, d, lb, ub, runge)

		pred_cheb = predict(cheb, points)
		pred_uni = predict(unic, points)

		cheb_m[:, i] = pred_cheb["preds"]
		uni_m[:, i] = pred_uni["preds"]
	end

	subplot(121)
	PyPlot.plot(points, runge.(points), label="True")
	for (i, d) in enumerate(deg)
		PyPlot.plot(points, cheb_m[:, i,], label=d)
	end
	title("Chebyshev nodes")
	legend()
	subplot(122)
	PyPlot.plot(points, runge.(points), label="True")
	for (i, d) in enumerate(deg)
		PyPlot.plot(points, uni_m[:, i,], label=d)
	end
	title("Uniform nodes")
	legend()
	PyPlot.savefig(joinpath(dirname(@__FILE__),"..","q5.png"))

end



function q6()

	# compare 2 knot vectors with runge's function

	PyPlot.savefig(joinpath(dirname(@__FILE__),"..","q6.png"))

end

function q7()
	f(x) = abs(x)^0.5
	precision = 63
	x = range(-1, 1, length=precision)
	myknots = vcat(range(-1,stop = -0.1,length = 5), 0, 0, 0, range(0.1,stop = 1,length =5))

	y = f.(x)
	bs = BSpline(13, 3, -1, 1)
	bs2 = BSpline(myknots, 3)
	B = Array(getBasis(collect(x), bs))
	B2 = Array(getBasis(collect(x), bs2))
	c = B \ y
	c2 = B2 \ y
	subplot(311)
	PyPlot.plot(x, y, label="True function")
	PyPlot.scatter(x, B * c, color="red", s=2, label="Approximation")
	title("Uniform knot vector")
	legend()
	subplot(312)
	PyPlot.plot(x, y, label="True function")
	PyPlot.scatter(x, B2 * c2, color="red", s=2, label="Approximation")
	title("Knot multiplicity x=0")
	legend()
	subplot(313)
	PyPlot.plot(x, y - B * c, color="blue", label="Uniform")
	PyPlot.plot(x, y - B2 * c2, color="red", label="Knot multiplicity")
	title("Error plots")
	legend()
	PyPlot.savefig(joinpath(dirname(@__FILE__),"..","q7.png"))
end


	# function to run all questions
function runall()
	@info("running all questions of HW-funcapprox:")
	q1(15)
	q2(3)
	q3(10)
	q4()
	q5()
	q6()
	q7()
end



end # module
