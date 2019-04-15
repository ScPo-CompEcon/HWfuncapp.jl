module HWfuncapp

using FastGaussQuadrature # to get chebyshevnodes

# you dont' have to use PyPlot. I did much of it in PyPlot, hence
# you will see me often qualify code with `PyPlot.plot` or so.
using PyPlot
import ApproXD: getBasis, BSpline
using Distributions
using ApproxFun
using Plots

export scalemap, unitmap, ChebyT, q1, q2, q3, runall

ChebyT(x,deg) = cos(acos(x)*deg)
unitmap(x,lb,ub) = 2 .* (x .- lb) ./ (ub .- lb) .- 1	#[a,b] -> [-1,1]
scalemap(x,lb,ub) = 0.5 .* (lb .+ ub) .+ 0.5 .* (ub .- lb) .* x     #from [-1,1] to [a,b]

function q1(n)
    f(x) = x .+ 2x.^2 - exp.(-x)
    lb = -3
    ub = 3
    n_new = 100
    grid = range(lb, length = n_new, stop = ub)
    chNodes = gausschebyshev(n)[1]

    phi = zeros(length(chNodes), n)
    for i in 1:length(chNodes)
        for j in 1:n
            phi[i,j] = ChebyT(chNodes[i], j-1)
        end
    end
    c = inv(phi)*f(scalemap(chNodes, lb, ub))

    true_values = f(grid)
    approx_matrix = zeros(length(grid),n)
    for i in 1:length(grid)
        for j in 1:n
            approx_matrix[i,j] = ChebyT(unitmap(grid,lb,ub)[i], j-1)
        end
    end
    approx_values = approx_matrix*c
    error = true_values - approx_values

    p = Plots.plot(layout = 2)
    Plots.plot!(p[1], grid, true_values, label = "f(x)")
    Plots.plot!(p[1], grid, approx_values, label = "approx f(x)", linestyle = :dot)
    Plots.plot!(p[2], grid, error, label = "error")
    Plots.savefig(p,joinpath(dirname(@__FILE__),"..","q1.png"))
    return Dict(:error=>maximum(abs,error))
end

function q2(b::Number)
    @assert b > 0
    f_true = Fun(x -> x .+ 2x.^2 - exp.(-x),-b..b)
    n = 100
    grid = range(-b, length = n, stop = b)
    true_values = f_true.(grid)

    cheby = Chebyshev(-b..b)
    diagonal = zeros(n,n)
    for i = 1:n
        diagonal[:,i] = Fun(cheby,[zeros(i-1);1]).(grid)
    end
    f_app = Fun(cheby,diagonal\true_values)
    approx_values = f_app.(grid)
    error = true_values - approx_values

    p = Plots.plot(layout = 2)
    Plots.plot!(p[1], grid,true_values, label = "f(x)")
    Plots.plot!(p[1], grid, approx_values, label = "approx f(x)", linestyle = :dot)
    Plots.plot!(p[2], grid, error, label = "error")
    Plots.savefig(p,joinpath(dirname(@__FILE__),"..","q2.png"))
end

function q3(b::Number)
    @assert b > 0
	x = Fun(identity, -b..b)
    f = sin(x^2)
    g = cos(x)
    h = f - g
    roots_h = roots(h)
	a = cumsum(h)
	integral = a(0) - a(-b)

    p = Plots.plot(h, label = "h(x)")
    Plots.scatter!(p , roots_h , h.(roots_h), label = "roots h(x)")
    Plots.savefig(p , joinpath(dirname(@__FILE__),"..","q3.png"))

    return (p, integral)
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
