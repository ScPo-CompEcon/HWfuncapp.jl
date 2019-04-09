module HWfuncapp

using FastGaussQuadrature # to get chebyshevnodes

# you dont' have to use PyPlot. I did much of it in PyPlot, hence
# you will see me often qualify code with `PyPlot.plot` or so.
using PyPlot
import ApproXD: getBasis, BSpline
using Distributions
using BasisMatrices
using ApproxFun
using Plots

export q1, q2, q3, q7

ChebyT(x,deg) = cos(acos(x)*deg)
unitmap(x,lb,ub) = 2 .* (x .- lb) ./ (ub .- lb) .- 1	#[a,b] -> [-1,1]
abmap(x,lb,ub) = 0.5 .* (ub .+ lb) .+ 0.5 .* (ub .- lb) .* x	# [-1,1] -> [a,b]
function chebpol(x, n)
	Φ = zeros(size(x,1),n) # chebyshev polynomials basis
	for i in 1:size(x,1)
		for j in 1:n
			Φ[i,j] = ChebyT(x[i],j-1)
			# Φ[i,j] = abmap(Φ[i,j], lb, ub)
		end
	end
	return Φ
end


function q1(n=15)
	lb = -3
	ub = 3
	x = range(lb, stop = ub ,length = n)
	# generate Chebyshev nodes
	S, y = gausschebyshev(n)
	# Scale Chebyshev nodes on -3 to 3 domain
	ϕ = abmap(S, lb, ub) # Chebyshev nodes on -3, 3
	# define function
	f(x) = x .+ 2x.^2 - exp.(-x)
	# Evaluate function at chebyshev nodes
	Y = f(ϕ)

	# Build Function grid
	Φ = chebpol(S, n)

	# Estimate parameters of approximation
	c = Φ\Y

	# Evaluate approximation
	Yhat = Φ*c

	# Test accuracy on larger sample
	n_new = range(lb, stop = ub ,length = 100)
	# Evaluate function
	Y_new = f(n_new)
	# Build Polynomial grid
	Φ_new = chebpol(unitmap(n_new,lb,ub), n)

	# Evaluate using fitted coefficients
	Yhat_new = Φ_new*c

	# Plot
	p = Any[]
	push!(p,Plots.plot(n_new, [Y_new, Yhat_new], label=["True value" "Approximation"], marker=([:none :diamond])))
	err = Y_new .- Yhat_new
	push!(p, Plots.plot(n_new, err, label="Approximation error"))
	Plots.plot(p...)
	# without using PyPlot, just erase the `PyPlot.` part
	Plots.savefig(Plots.plot(p...), joinpath(dirname(@__FILE__),"..","q1.png"))
	return Dict("error"=>maximum(abs,err))
end


function q2(b=4)
	@assert b > 0
	n = 15
	# use ApproxFun.jl to do the same:
	deg = n - 1
	ub = b
	lb = -b
	S = Chebyshev(lb..ub)
	p = range(lb,stop=ub,length=n)  # a non-default grid
	# define function
	f(x) = x .+ 2x.^2 - exp.(-x)
	# evaluate function
	v = f(p)           # values at the non-default grid
	V = Array{Float64}(undef,n,deg + 1) # Create a Vandermonde matrix by evaluating the basis at the grid

	for k = 1:deg+1
	    V[:,k] = Fun(S,[zeros(k-1);1]).(p)
	end
	V
	g = Fun(S,V\v);
	@show g(1.1)
	@show f(1.1)

	n_new = range(lb, stop = ub, length = 100)
	p = Any[]
	push!(p,Plots.plot(n_new, [f(n_new), g.(n_new)], label=["True value" "Approximation"], marker=([:none :diamond])))
	err = f(n_new) .- g.(n_new)
	push!(p, Plots.plot(n_new, err, label="Approximation error"))
	fig = Plots.plot(p...)

	Plots.savefig(fig,joinpath(dirname(@__FILE__),"..","q2.png"))
end

function q3(b::Number)

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
