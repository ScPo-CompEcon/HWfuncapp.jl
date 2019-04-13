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

function q1(n)
	f(x) = x .+ 2x.^2 - exp.(-x)
# Generate Chebychev nodes
	deg = n:-1:1
	z = 3*cos.((2 .* deg .- 1).*pi ./ 2n)

	# Compute function values at interpolation nodes
	y = zeros(n)
	for i in 1:n
       y[i] = f(z[i])
    end

	# Compute basis matrix at the nodes
	P = zeros(n, n)
	for i in 0:n-1
	       for j in 1:n
	       P[j,i+1] = ChebyT.(unitmap(z, -3, 3)[j], i)
	       end
	end
	# Obtain coefficients (solve the system)
	c = inv(P) * y

	# Recompute matrix for new set of points
	# Evaluate Φ(x'): i.e. evaluate all 15 Chebyshev polynomials at the new x, and multiply them by c
	n_new = 100
	b = unitmap(range(-3, stop = 3, length = n_new), -3, 3)
	P_b = zeros(n_new, n)
	for i in 0:n-1
	       for j in 1:n_new
	       P_b[j,i+1] = ChebyT(b[j], i)
	       end
	end

	# Multiply this new matrix by c to obtain approximation of f
		y_new = f.(range(-3, stop = 3, length = n_new))
		y_approx = (P_b) * c

		approx_error = y_new .- y_approx
		plot(range(-3, stop = 3, length = n_new), [y_new approx_error], label = "Actuals", layout = 2)
		scatter!(range(-3, stop = 3, length = n_new), y_approx, label = "Approximation", marker = (1, :circle))

	# without using PyPlot, just erase the `PyPlot.` part
	savefig(joinpath(dirname(@__FILE__),"..","q1.png"))
	return Dict(:error=>maximum(abs,approx_error))
end

function q2(b::Number)
	@assert b > 0
	# use ApproxFun.jl to do the same:

	Plots.savefig(p,joinpath(dirname(@__FILE__),"..","q2.png"))
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
