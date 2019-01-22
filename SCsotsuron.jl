using Base.Test
import Base.digits
using Plots; plotly()
using NLsolve # for solving the nonliniar equations
using SincFun, DEQuadrature # for running DEQuadrature

include("Getparameter.jl")
include("PlotFunctions.jl")
include("SCsupport.jl")
include("Trapezoidal.jl")

# ---------------------------------------------------------------------------------------------------- #
### settings ###

## setting precision of BigFloat
digits(n::Int) = setprecision(round(Int,ceil(n*log2(10))))
digits(300)

## poles and the integral function
## Example 4.1 from Slevinsky et al.(2015)

#=
exid = 10
z = [complex(-0.5,1.0);complex(0.5,0.5)]
NpSL = 2
dlt = zto_atanh(z)[1]
ept = zto_atanh(z)[2]
Np = length(dlt)
f = x-> exp.(1./abs2.(x-z[1])).*log.(1-x)./abs2.(x-z[2])./sqrt.(1+x)
T = 0.5*log(big(2.0))
β_factorDE = big(0.50)
β_factorSL = big(0.50)
β_factorNN = 1/sqrt(big(2.0))
d_DE = minimum(imag.(asinh.((2/big(π)).*atanh.(z))))
d_SL = big(π)/2
d_NN = big(π)/2
d_N2 = big(π)/2
Itrue = "-2.04645081160694748690442050179886173463698400851312978159495108281833923937599915411239665314732909414150138401559892277212580897054853305044117124435121918549124604645789190248721167573993297627891884242759929922470299331035480062586543433278436544374899377193299256405496172603459957078880800178958528771642962893695508317394167405216487683527096100475182610564645687173201068270992551244010304597e+00"
=#

## example(番号1) : Example4.2 from Slevinsky et al.(2015)

#=
exid = 20
z = [complex(big(-2.0),1.0);
     complex(-1.0,0.5);
     complex(1.0,0.25);
     complex(2.0,1.0)]
NpSL = 4
Np = 4
dlt = zto_asinh(z)[1]
ept = zto_asinh(z)[2]
f = x-> exp.(10./abs2.(x-z[1])).*cos.(10./abs2.(x-z[2]))./abs2.(x-z[3])./abs.(x-z[4])
T = big(0.0)
β_factorDE = big(1.0)
β_factorSL = big(1.0)
β_factorNN = big(1.0)
d_DE = minimum(imag.(asinh.((2/big(π)).*asinh.(z))))
d_SL = big(π)/2
d_NN = big(π)/2
d_N2 = big(π)/2
Itrue = "1.50133619876062770101030470326173553208854739646240081225845195322624377330867094144192464692844931351956249439446687259575908821074467954700173593866783934590906808530844595015500541680530415264580194453133658724574536164866695611699188173914168980057456962540582381978462005370937587010012338985572909021270993977527355156445163416743376725727256310809095229447802636632237845095430431700882602931e+01"#
=#

# ひどい例
#=
exid = 300
z = [complex(big(-2.0),0.25);complex(-1.5,0.25);complex(-1.0,0.5);complex(0.0,0.5);
complex(0.5,0.5);complex(1.5,0.25);complex(2.0,0.25)]
dlt = zto_asinh(z)[1]
ept = zto_asinh(z)[2]
Np = length(dlt)
f = x -> exp.(1./abs2.(x-z[3])).*exp.(1./abs2.(x-z[5])).*cos.(10./abs2.(x-z[1])).*cos.(5./abs2.(x-z[7]))./abs.(x-z[2])./abs.(x-z[4])./abs.(x-z[6])
Itrue = "-2.55424501539945020894237800465334255132690317209188447333433600496112290647081743581978940284053825625408339404217315475470816794345781139993960918433496419279049768486132511050089232890815668766452007"
β_factor = 1
T = 0
=#

# ひどい例その2
#=
exid = 400
z = [complex(big(0.2),0.20);complex(0.6,0.1);
complex(1.0,0.5);complex(1.5,0.25);
complex(2.0,0.5);complex(2.5,0.25);complex(3.0,0.25)]
zzz = [complex(big(0.2),0.20);complex(3.0,0.25)]
#dlt = zto_log(z)[1]
#ept = zto_log(z)[2]
dlt = zto_log(zzz)[1]
ept = zto_log(zzz)[2]
Np = length(dlt)
f = x -> cos.(2./abs2.(x-z[1])).*exp.(1./abs2.(x-z[3])).*exp.(1./abs2.(x-z[5])).*cos.(1./abs2.(x-z[7]))./sqrt.(x)./abs.(x-z[2])./abs.(x-z[4])./abs.(x-z[6])
Itrue = "-1.039273975787783078395741438068104103249411393444261027746957420793688463925546037330823845841155525229547587450307910969437396515722364442824828420423534359648919639293555052850123094783735102170087574e+01"
T = -log(big(2.0))
β_factorSL = 0.0
β_factorDE = big(0.25)
β_factorNN = big(0.50)
=#

#　Slevinsky がギリギリいける
#=
exid = 500
z = [complex(big(0.5),0.1);complex(0.8,1.0);complex(1.0,0.5);complex(1.5,0.3);complex(2.0,0.6)]
f = x -> cos.(2./abs2.(x-z[1])).*exp.(1./abs2.(x-z[3])).*exp.(1./abs2.(x-z[5]))./sqrt.(x)./abs.(x-z[2])./abs2.(x-z[4])
Itrue = "-5.457407670751856530484750347002849615931930095181485083457494004640400248846376222846032549075905723750063650054306682999136362687471636222480400261690929903294253382948508434017142728587658884021624639e+01"
dlt = zto_log(z)[1]
ept = zto_log(z)[2]
NpSL = 5
Np = length(dlt)
T = -log(big(2.0))
β_factorDE = big(0.25)
β_factorSL = big(0.25)
β_factorNN = big(0.50)
d_DE = minimum(imag.(asinh.((2/big(π)).*log.(z))))
# d_DE = big(π)/2
d_SL = big(π)/2
d_NN = big(π)/2
d_N2 = big(π)/2
=#

# Slevinsky の単射性がなくなる
#=
exid = 600
z = [complex(big(0.3),0.2);complex(0.5,0.6);complex(0.8,0.5);complex(1.2,0.3)]
f = x -> exp.(0.02./abs.(x-z[1]).^3).*exp.(0.05./abs.(x-z[4]).^3)./sqrt.(x)./abs.(x-z[2])./abs2.(x-z[3])
Itrue = "3.069293691227855355300418056148943017596912116830751277390412174563778418345716563374389489803603308701717053926393928160305389745320089440166780387281514261598188082078924668997650240748322997832870553e+01"
NpSL = 4
dlt = zto_log(z)[1]
ept = zto_log(z)[2]
Np = 4
T = -log(big(2.0))
β_factorDE = big(0.25)
β_factorSL = big(0.25)
β_factorNN = big(0.50)
d_DE = minimum(imag.(asinh.((2/big(π)).*log.(z))))
d_SL = big(π)/2
d_NN = big(π)/2
d_N2 = big(π)/2
=#

#あまりにもひどい例

exid = 1000
z = [complex(big(1.0),0.1);complex(2.0,0.5);complex(3.0,0.3);complex(4.0,0.5);
     complex(5.0,0.2);complex(6.0,0.5);complex(7.0,0.1)]
f = x -> @. cos(5/abs2(x-z[1]))*cos(10/abs2(x-z[7]))*
            exp(0.8/abs2(x-z[2]))*exp(0.2/abs2(x-z[3]))*
            exp(0.5/abs2(x-z[4]))*exp(0.1/abs2(x-z[5]))*exp(0.5/abs2(x-z[6]))*
            exp(-0.2x)/sqrt(x)
Itrue = "-3.451882594217500854747246631701525068156038028658631544293802399002537897120929315856041626221574912e-01"
dlt = zto_logexp(z)[1]
ept = zto_logexp(z)[2]
Np = length(dlt)
T = - log(big(2.5))/2
β_factorDE = big(0.10)
β_factorNN = sqrt(big(0.1))/2
d_DE = minimum(imag.(asinh.((2/big(π)).*log.(exp.(z)-1))))
d_NN = big(π)/2
β_factorSL = 0.0
d_SL = 0.0


# ---------------------------------------------------------------------------------------------------- #
### the DE transformation ###
hDE  = t -> (big(π)./2).*sinh.(t)
dhDE = t -> (big(π)./2).*cosh.(t)
FDEsinh  = t -> f(sinh.(hDE(t))).*cosh.(hDE(t)).*dhDE(t)
FDEtanh = t -> f(tanh.(hDE(t))).*dhDE(t)./(cosh.(hDE(t))).^2
FDEexp = t -> f(exp.(hDE(t))).*exp.(hDE(t)).*dhDE(t)
FDElogexp = t -> f(log.(exp.(hDE(t))+1)).*exp.(hDE(t)).*dhDE(t)./(exp.(hDE(t))+1)
β_DE  = (big(π)/2)*β_factorDE

# ---------------------------------------------------------------------------------------------------- #
### the method of Slevinsky et al.(2015) ###

print("getting parameters: Slevinsky\n")
if exid == 10
    h = DEMapValues(z;domain=Finite(-1.0,1.0,-0.5,0.0,0.0,1.0))
    @time h  = DEMapValues(z;domain=Finite(-1.0,1.0,-0.5,0.0,0.0,1.0))
end
if exid == 20
    h  = DEMapValues(z;domain=Infinite1(BigFloat),Hint=10,obj_scaling_factor=-1.0)
    @time h  = DEMapValues(z;domain=Infinite1(BigFloat),Hint=10,obj_scaling_factor=-1.0)
end
if exid == 300 || exid == 400 || exid == 1000
    # 制御するのがめんどいので無意味なやつを入れておく
    h = DEMapValues([big(1.0)+2.0im];domain=Infinite1(BigFloat),Hint=10,obj_scaling_factor=-1.0)
end
if exid == 500
    h = DEMapValues(z;domain=SemiInfinite2(BigFloat),Hint=10000,obj_scaling_factor=-1.0)
    @time h = DEMapValues(z;domain=SemiInfinite2(BigFloat),Hint=10000,obj_scaling_factor=-1.0)
end
if exid == 600
    h = DEMapValues(z;domain=SemiInfinite2(BigFloat),Hint=10000,obj_scaling_factor=-0.5)
    #@time h = DEMapValues(z;domain=SemiInfinite2(BigFloat),Hint=10000,obj_scaling_factor=-0.5)
end
u0 = h.u0
u  = h.u
function hSL(t)
    f = u0.*sinh.(t)
    for i = 1:NpSL
        f = f + u[i].*t.^(i-1)
    end
    return f
end
function dhSL(t)
    df = u0.*cosh.(t)
    for i=2:NpSL
        df = df + (i-1).*u[i].*t.^(i-2)
    end
    return df
end
FSLsinh  = t -> f(sinh.(hSL(t))).*cosh.(hSL(t)).*dhSL(t)
FSLtanh = t -> f(tanh.(hSL(t))).*dhSL(t)./(cosh.(hSL(t))).^2
FSLexp = t -> f(exp.(hSL(t))).*exp.(hSL(t)).*dhSL(t)
β_SL  = u0*β_factorSL

# ---------------------------------------------------------------------------------------------------- #
### the proposed method ###

## getting paramerers of the proposed formula
D = zeros(BigFloat,Np-1)
for i = 1:Np-1
    D[i] = (dlt[i+1]-dlt[i])/big(π)
end
print("getting parameters: OptSC\n")
Cab = xtoCab(Getparameter().zero)
@time Cab = xtoCab(Getparameter().zero)
#PlotFimag(Cab,ept) # plotting the state of convergence
C = Cabsprit(Cab)[1]
a = Cabsprit(Cab)[2]
b = Cabsprit(Cab)[3]

## the proposed transformation formula
function hNN(t)
    f = C.*sinh.(t-T) + dlt[1]
    for i = 1:Np-1
        f = f + 2.*D[i].*atan.(exp.(t-b[i]))
    end
    return f
end
function dhNN(t)
    df = C.*cosh.(t-T)
    for i = 1:Np-1
        df = df + D[i]./cosh.(t-b[i])
    end
    return df
end
function hN2(t)
    f = C.*sinh.(t-T) + (dlt[1]+dlt[Np])./2
    for i = 1:Np-1
        f = f + D[i].*(big(π)./2.*tanh.((t-b[i]).*2./big(π)))
    end
    return f
end
function dhN2(t)
    df = C.*cosh.(t-T)
    for i = 1:Np-1
        df = df + D[i]./cosh.((t-b[i]).*2./big(π)).^2
    end
    return df
end

FNNsinh = t -> f(sinh.(hNN(t))).*cosh.(hNN(t)).*dhNN(t) #example(番号1)
FNNtanh = t -> f(tanh.(hNN(t))).*dhNN(t)./(cosh.(hNN(t))).^2
FNNexp = t -> f(exp.(hNN(t))).*exp.(hNN(t)).*dhNN(t)
FNNlogexp = t -> f(log.(exp.(hNN(t))+1)).*exp.(hNN(t)).*dhNN(t)./(exp.(hNN(t))+1)
FN2sinh = t -> f(sinh.(hN2(t))).*cosh.(hN2(t)).*dhN2(t) #example(番号1)
FN2tanh = t -> f(tanh.(hN2(t))).*dhN2(t)./(cosh.(hN2(t))).^2
FN2exp = t -> f(exp.(hN2(t))).*exp.(hN2(t)).*dhN2(t)
FN2logexp = t -> f(log.(exp.(hN2(t))+1)).*exp.(hN2(t)).*dhN2(t)./(exp.(hN2(t))+1)
β_NN  = C*β_factorNN

# ---------------------------------------------------------------------------------------------------- #
### comparing the methods ###

## conformal mapping
# PlotConformal(hN2, -15.0:0.01:15, 5)
# PlotConformal(hSL, -18.0:0.1:15, 5)
# PlotConformal_b(hNN,5,a,b) # conformal mapping of hNN

## plotting transformed integral functions

p1 = plot()
xx = -10.0:0.01:10
if exid == 10
    plot!(p1, xx, FSLtanh(xx), color="blue", label="SO", line=:dashdot)
    plot!(p1, xx, FDEtanh(xx), color="green", label="DE")
    plot!(p1, xx, FN2tanh(xx), color="yellow", label="New2", w=3)
    plot!(p1, xx, FNNtanh(xx), color="red", label="New", w=3, line=:dot)
end
if exid == 20
    # plot!(p1, xx, f(xx), color="black", label="Original")
    plot!(p1, xx, FSLsinh(xx), color="blue", label="SO", line=:dashdot)
    plot!(p1, xx, FDEsinh(xx), color="green", label="DE")
    plot!(p1, xx, FN2sinh(xx), color="yellow", label="New2", w=3)
    plot!(p1, xx, FNNsinh(xx), color="red", label="New", w=3, line=:dot)
end
if exid == 300
    plot!(p1, xx, f(xx), color="black", label="Original")
    plot!(p1, xx, FDEsinh(xx), color="green", label="DE")
    plot!(p1, xx, FNNsinh(xx), color="red", label="New")
    plot!(p1, xx, FN2sinh(xx), color="yellow", label="New2")
end
if exid == 400
    plot!(p1, xx, FDEexp(xx), color="green", label="DE")
    plot!(p1, xx, FNNexp(xx), color="red", label="New")
    plot!(p1, xx, FN2exp(xx), color="yellow", label="New2")
end
if exid == 500 || exid == 600
    plot!(p1, xx, FSLexp(xx), color="blue", label="SO", line=:dashdot)
    plot!(p1, xx, FDEexp(xx), color="green", label="DE")
    plot!(p1, xx, FN2exp(xx), color="yellow", label="New2", w=3)
    plot!(p1, xx, FNNexp(xx), color="red", label="New", w=3, line=:dot)
end
if exid == 1000
    plot!(p1, xx, FDElogexp(xx), color="green", label="DE")
    plot!(p1, xx, FN2logexp(xx), color="yellow", label="New2", w=3)
    plot!(p1, xx, FNNlogexp(xx), color="red", label="New", w=3, line=:dot)
end
display(p1)


## numerical integralation
## setting the order n
if exid == 10
    nmax = 10
    nidx = [30*i for i=1:nmax]
end
if exid == 20
    nmax = 20
    nidx = [30*i for i=1:nmax]
end
if exid == 300
    nmax = 10
    nidx = [2^i for i=1:nmax]
end
if exid == 400
    nmax = 20
    nidx = [30*i for i=1:nmax]
end
if exid == 500
    nmax = 10
    nidx = [100*i for i=1:nmax]
end
if exid == 600
    nmax = 20
    nidx = [30*i for i=1:nmax]
end
if exid == 1000
    nmax = 10
    nidx = [500*i for i=1:nmax]
end
IDE = zeros(BigFloat,nmax)
INN = zeros(BigFloat,nmax)
ISL = zeros(BigFloat,nmax)
IN2 = zeros(BigFloat,nmax)
if exid == 10
    nIarray!([0.0],FDEtanh,[2],β_DE,d_DE)
    print("integral: FDE\n")
    @time nIarray!(IDE,FDEtanh,nidx,β_DE,d_DE)
    print("integral: FSL\n")
    @time nIarray!(ISL,FSLtanh,nidx,β_SL,d_SL)
    print("integral: FNN\n")
    @time nIarray!(INN,FNNtanh,nidx,β_NN,d_NN)
    print("integral: FN2\n")
    @time nIarray!(IN2,FN2tanh,nidx,β_NN,d_NN)
end
if exid == 20
    nIarray!([0.0],FDEtanh,[2],β_DE,d_DE)
    print("integral: FDE\n")
    @time nIarray!(IDE,FDEsinh,nidx,β_DE,d_DE)
    print("integral: FSL\n")
    @time nIarray!(ISL,FSLsinh,nidx,β_SL,d_SL)
    print("integral: FNN\n")
    @time nIarray!(INN,FNNsinh,nidx,β_NN,d_NN)
    print("integral: FN2\n")
    @time nIarray!(IN2,FN2sinh,nidx,β_NN,d_N2)
end
if exid == 300
    nIarray!([0.0],FDEtanh,[2],β_DE,d_DE)
    print("integral: FDE\n")
    @time nIarray!(IDE,FDEsinh,nidx,β_DE)
    print("integral: FNN\n")
    @time nIarray!(INN,FNNsinh,nidx,β_NN)
    print("integral: FN2\n")
    @time nIarray!(IN2,FN2sinh,nidx,β_NN)
end
if exid == 400
    nIarray!([0.0],FDEtanh,[2],β_DE,d_DE)
    print("integral: FDE\n")
    @time nIarray!(IDE,FDEexp,nidx,β_DE)
    print("integral: FNN\n")
    @time nIarray!(INN,FNNexp,nidx,β_NN)
    print("integral: FN2\n")
    @time nIarray!(IN2,FN2exp,nidx,β_NN)
end
if exid == 500 || exid == 600
    nIarray!([0.0],FDEexp,[2],β_DE,d_DE)
    print("integral: FDE\n")
    @time nIarray!(IDE,FDEexp,nidx,β_DE,d_DE)
    print("integral: FSL\n")
    @time nIarray!(ISL,FSLexp,nidx,β_SL,d_SL)
    print("integral: FNN\n")
    @time nIarray!(INN,FNNexp,nidx,β_NN,d_NN)
    print("integral: FN2\n")
    @time nIarray!(IN2,FN2exp,nidx,β_NN,d_N2)
end
if exid == 1000
    nIarray!([0.0],FDElogexp,[2],β_DE,d_DE)
    print("integral: FDE\n")
    @time nIarray!(IDE,FDElogexp,nidx,β_DE,d_DE)
    print("integral: FNN\n")
    @time nIarray!(INN,FNNlogexp,nidx,β_NN,d_NN)
    print("integral: FN2\n")
    @time nIarray!(IN2,FN2logexp,nidx,β_NN,d_NN)
end

## comparing results of numerical integration
# p2 = plot()
# plot!(p2, nidx, IDE, color="green", label="DE")
# plot!(p2, nidx, ISL, color="blue", label="SL")
# plot!(p2, nidx, INN, color="red", label="New")
# hline!([parse(BigFloat,Itrue)], color="black",label="")
# display(p2)

# EEmax = n-> @. -π^2*(2n+1)/(2*log(π^2*(2n+1)/2/β_SL)*log(10))
## comparing the relative error
myerr = x -> log10.(abs.((x-parse(BigFloat,Itrue))./parse(BigFloat,Itrue)))
p3 = plot()
if exid == 10 || exid == 20 || exid == 500 || exid == 600
    plot!(p3, nidx, myerr(ISL), color="blue", label="SO", marker=(3,:rect), line=:dashdot)
end
plot!(p3, nidx, myerr(IDE), color="green", label="DE", marker=(3,:circle))
plot!(p3, nidx, myerr(IN2), color="yellow", label="New2", marker=(4,:dtriangle), w=3)
plot!(p3, nidx, myerr(INN), color="red", label="New", marker=(4,:utriangle), w=3, line=:dot)
xlabel!(p3, "order n")
ylabel!(p3, "log10|RelativeError|")
display(p3)

# display(plot(dlt,ept,seriestype=:scatter))

# Float16.(imag.(asinh.((2./pi).*asinh.(z))))
