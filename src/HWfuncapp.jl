module HWfuncapp

using FastGaussQuadrature # to get chebyshevnodes

# you dont' have to use PyPlot. I did much of it in PyPlot, hence
# you will see me often qualify code with `PyPlot.plot` or so.
using PyPlot
import ApproXD: getBasis, BSpline
using Distributions
using ApproxFun
using Plots

ChebyT(x,deg) = cos(acos(x)*deg)
unitmap(x,lb,ub) = 2 .* (x .- lb) ./ (ub .- lb) .- 1	#[a,b] -> [-1,1]

function q1(n = 15)
	f(x) = x .+ 2x.^2 - exp.(-x)
	deg = n-1
	lb = -3
	ub = 3
	node = gausschebyshev(n)[1]
	chebnode = 0.5 .* (lb .+ ub) .+ 0.5 .* (ub .- lb) .* node
	y = map(f, chebnode)
	basismatrix = [cos((n - i + 0.5) * (j - 1)* π / n) for i = 1:n, j = 1:n]
	c = basismatrix \ y

	n_new = 100
	x_2 = range(lb, stop = ub, length = n_new)
	base2 = [ChebyT(unitmap(x_2[i], lb, ub), j) for i in 1:n_new, j = 0:deg]
	y2 = base2 * c
	error = y2 - f(collect(-3:6/99:3))
	Plots.plot(x_2, [f(collect(-3:6/99:3)) error], labels = ["true_val" "error"], title = "Q1", layout = 2)
	scatter!(x_2, y2, labels = ["approx"])

	# without using PyPlot, just erase the `PyPlot.` part
	Plots.savefig(joinpath(dirname(@__FILE__),"..","q1.png"))
	return Dict(:error=>maximum(abs,error))
end

function q2(b::Number)
	@assert b > 0
	# use ApproxFun.jl to do the same:
	b = 4
	n_new = 100
	f(x) = x .+ 2x.^2 - exp.(-x)
	f2 = Fun(f,-b..b)
	x = range(-b, b, length = n_new)
	y2 = f2.(x)
	y = f.(x)
	error = y - y2
	Plots.plot(x, [y error], labels = ["true_val" "error"], title = "Q2", layout = 2)
	scatter!(x, y2, labels = ["approx"])

	Plots.savefig(joinpath(dirname(@__FILE__),"..","q2.png"))
end
using LinearAlgebra #for cumsum
function q3(b::Number)
	b = 10
	x = Fun(identity, -b..b)
	f = sin(x^2)
	g = cos(x)
	h = f - g
	r = roots(h)
	p = Plots.plot(h, label = "h", title = "Q3")
	scatter!(r,h.(r),labels="roots h")
	Plots.savefig(joinpath(dirname(@__FILE__),"..","q3.png"))
	Plots.savefig(joinpath(p,dirname(@__FILE__),"..","q3.png"))
	k1 = cumsum(h)
	k = k1 + h(-b)
	integral = norm(h-k)
	# p is your plot
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

function q7()


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
