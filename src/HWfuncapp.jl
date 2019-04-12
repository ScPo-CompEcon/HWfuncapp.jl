module HWfuncapp

# you dont' have to use PyPlot. I did much of it in PyPlot, hence
# you will see me often qualify code with `PyPlot.plot` or so.
import ApproXD: getBasis, BSpline
using Distributions
using ApproxFun
using Plots
using BasisMatrices
using FastGaussQuadrature
using LinearAlgebra, SpecialFunctions
export ChebyT,Chebypolynomial, unitmap, q1, q2, q3, q7

ChebyT(x,deg) = cos(acos(x)*deg)
unitmap(x,lb,ub) = 2 .* (x .- lb) ./ (ub .- lb) .- 1	#[a,b] -> [-1,1]
reversemap(x,lb,ub) = 0.5 .* (x .+ 1) .* (ub .- lb) .+lb  #[-1,1] ->[a,b]
function Chebypolynomial(x,n)
	M = zeros(length(x),n)
	for i in 1:length(x)
		for j in 1:n
			M[i,j] = ChebyT(x[i],j-1) #the first column is ones
		end
	end
	return M
end

function q1(n=15)
	lb = -3
	ub = 3
	x = range(lb,stop=ub,length=n)
	S = gausschebyshev(n)[1]
	f(x) = x .+ 2x.^2 - exp.(-x)
	# evaluate function at Chebyshev nodes
	Y = f(reversemap(S,lb,ub))
    # get the chebyshev polynomial
    M = Chebypolynomial(S,n)
	# get the coeffcients
	c = inv(M)*Y
	# test 1
	n_new = 100
	x_new = range(-3,stop=3,length=n_new)
	#S_new, y_new = gausschebyshev(n_new)
	#Y_new = f(reversemap(S_new,lb,ub))
	#M_new = Chebypolynomial(S_new,n)
	Y_new = f(x_new)
	M_new = Chebypolynomial(unitmap(x_new,lb,ub),n)
	Yhat_new = M_new*c
	err = Y_new .- Yhat_new
	p = Plots.plot(layout = 2, dpi = 400)
	Plots.plot!(p[1],x_new,Y_new,label = "Y",lw = 1,linecolor = "black")
	Plots.plot!(p[1],x_new,Yhat_new, label = "Yhat", lw = 2,linestyle = :dot, linecolor = "red")
	Plots.plot!(p[2],x_new,err, label = "error", lw = 1, linecolor = "green")
	# without using PyPlot, just erase the `PyPlot.` part
	Plots.savefig(p,joinpath(dirname(@__FILE__),"..","q1.png"))
	return Dict(:error=>maximum(abs,err))
end

function q2(b::Number=4)
	@assert b > 0
	# use ApproxFun.jl to do the same:
	n = 15
	S = Chebyshev(-b..b)
	p = range(-b,stop=b,length=n)
	f = x->x .+ 2x.^2 - exp.(-x)
	Y = f.(p)
	V = zeros(n,n)
	for i in 1:n
		V[:,i] = Fun(S,[zeros(i-1);1]).(p)
	end
	funcapp = Fun(S,V\Y)
	Yhat = funcapp.(p)
	err = Yhat - Y
	q = Plots.plot(layout = 2, dpi = 400)
	Plots.plot!(q[1],p,Y,label = "Y",lw = 1,linecolor = "black")
	Plots.plot!(q[1],p,Yhat, label = "Yhat", lw = 1,linestyle = :dot, linecolor = "red")
	Plots.plot!(q[2],p,err, label = "error", lw = 2, linecolor = "green")
	Plots.savefig(q,joinpath(dirname(@__FILE__),"..","q2.png"))
	return Dict(:error=>maximum(abs,err))
end

function q3(b::Number=4)
    x = Fun(identity,-b..b)
	f = sin(x^2)
	g = cos(x)
	h = f-g
	rh = roots(h)
	rp = roots(h')
	# p is your plot
	p = Plots.plot(h, label = "h(x)", dpi = 400, lw = 1, linecolor = "black")
	Plots.scatter!(p,rh,h.(rh), label = "roots of h(x)", markercolor = "green", markersize = 5)
    Plots.scatter!(p,rp,h.(rp), label = "MinMax of h(x)", markercolor = "red", markersize = 5)
	Plots.savefig(p,joinpath(dirname(@__FILE__),"..","q3.png"))
	v1 = cumsum(h)
	v = v1 + f(-b)
	integral = norm(h-v)
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
