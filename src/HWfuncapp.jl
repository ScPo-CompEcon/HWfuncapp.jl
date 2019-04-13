module HWfuncapp

using FastGaussQuadrature # to get chebyshevnodes

# you dont' have to use PyPlot. I did much of it in PyPlot, hence
# you will see me often qualify code with `PyPlot.plot` or so.
using PyPlot  
import ApproXD: getBasis, BSpline
using Distributions
using ApproxFun
using Plots
export TT, ChebyT, q1


ChebyT(x,deg) = cos(acos(x)*deg)
unitmap(x,lb,ub) = 2 .* (x .- lb) ./ (ub .- lb) .- 1	#[a,b] -> [-1,1]

n=15
map(x,lb,ub) = 0.5 .* (lb .+ ub) .+ 0.5 .* (ub .- lb) .* x

function Cheby(x,deg)
len=length(x)
psi = ones(len,deg)
for i in 1:len
for j in 1:deg
psi[i,j] = ChebyT(x[i],j-1) 
end
end
return psi
end



function q1(n=15)
lb, ub, inter, Nod = -3, 3, 100, gausschebyshev(n)[1]#Values for evaluation
x = range(lb,stop = ub,length = n)#
psi = map(Nod, lb, ub) #Evaluation
    
f(x) = x .+ 2x.^2 - exp.(-x)#Function to evaluate
PSI=inv(Cheby(Nod,n))
ev=f(psi)
ResEv=PSI*ev
    
testX=range(lb, stop=ub, length=inter)
testY=f(testX)
   
testMap=map(testX, lb, ub)
testPsi=Cheby(unitmap(testMap), n)

ResEv=testPsi*PSI
    
testError=ev-ResEv
 
#PyPlot.plot()   
p = Plots.plot(layout = 2, dpi = 400)
Plots.plot!(p[1],testX,testY,label = "Test Y",lw = 1,linecolor = "black")
Plots.plot!(p[1],testX,ResEv, label = "Result Evaluation", lw = 2,linestyle = :dot, linecolor = "pink")
Plots.plot!(p[2],testX,testError, label = "Test Error", lw = 1, linecolor = "blue")
    
# without using PyPlot, just erase the `PyPlot.` part
PyPlot.savefig(joinpath(dirname(@__FILE__),"..","q1.png"))
return Dict(:error=>maximum(abs,err))
end


end # module



