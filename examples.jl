# "psid" corresponds to the function ψ.
# psid = 1: ψ = tanh(.)
#        2: ψ = sinh(.)
#        3: ψ = exp(.)
#        4: ψ = log(exp(.)+1)
#
# "singends" is the degree of endpoint singularities.
# if psid == 1
#     singends = [p; q] when f(x) = O((1-x)^p) (x →  1-0)
#                                 = O((1+x)^q) (x → -1+0)
# if psid == 2
#     singends = [r; s] when f(x) = O(|x|^r) (x → +∞)
#                                 = O(|x|^s) (x → -∞)
# if psid == 3
#     singends = [r; q] when f(x) = O(x^r) (x → +∞)
#                                 = O(x^q) (x → +0)
# if psid == 4
#     singends = [v; q] when f(x) = O(exp(-(v-ε)x)) (x → +∞)
#                                 = O(x^q) (x → +0)

using Base.Test
import Base.digits

### This program uses NLsolve for solving nonliniar equations.
using NLsolve

#################################################################
### setting precision of BigFloat
digits(n::Int) = setprecision(round(Int,ceil(n*log2(10))))
digits(300) # significant digits
#################################################################

include("Getparameters.jl")
include("Stilde.jl")
include("GetCab.jl")
include("HNew.jl")
include("Trapezoidal.jl")

exid = 4
print("Example: ", exid, "\n")

if exid == 1
    psid = 1
    S = [complex(big(-0.5),1.0); complex(0.5,0.5)]
    f = x-> exp.(1./abs2.(x-S[1])).*log.(1-x)./abs2.(x-S[2])./sqrt.(1+x)
    singends = [big(0.0); -0.5]
    Itrue = "-2.04645081160694748690442050179886173463698400851312978159495108281833923937599915411239665314732909414150138401559892277212580897054853305044117124435121918549124604645789190248721167573993297627891884242759929922470299331035480062586543433278436544374899377193299256405496172603459957078880800178958528771642962893695508317394167405216487683527096100475182610564645687173201068270992551244010304597e+00"
elseif exid == 2
    psid = 2
    S = [complex(big(-2.0),1.0);complex(-1.0,0.5);complex(1.0,0.25);complex(2.0,1.0)]
    f = x-> exp.(10./abs2.(x-S[1])).*cos.(10./abs2.(x-S[2]))./abs2.(x-S[3])./abs.(x-S[4])
    singends = [-big(3.0);-3.0]
    Itrue = "1.50133619876062770101030470326173553208854739646240081225845195322624377330867094144192464692844931351956249439446687259575908821074467954700173593866783934590906808530844595015500541680530415264580194453133658724574536164866695611699188173914168980057456962540582381978462005370937587010012338985572909021270993977527355156445163416743376725727256310809095229447802636632237845095430431700882602931e+01"
elseif exid == 3
    psid = 3
    S = [complex(big(0.3),0.2);complex(0.5,0.6);complex(0.8,0.5);complex(1.2,0.3)]
    f = x -> exp.(0.02./abs.(x-S[1]).^3).*exp.(0.05./abs.(x-S[4]).^3)./sqrt.(x)./abs.(x-S[2])./abs2.(x-S[3])
    singends = [-big(3.0);-0.5]
    Itrue = "3.069293691227855355300418056148943017596912116830751277390412174563778418345716563374389489803603308701717053926393928160305389745320089440166780387281514261598188082078924668997650240748322997832870553e+01"
elseif exid == 4
    psid = 4
    S = [complex(big(1.0),0.1);complex(2.0,0.5);complex(3.0,0.3);complex(4.0,0.5);complex(5.0,0.2);complex(6.0,0.5);complex(7.0,0.1)]
    f = x -> @. cos(5/abs2(x-S[1]))*cos(10/abs2(x-S[7]))*exp(0.8/abs2(x-S[2]))*exp(0.2/abs2(x-S[3]))*exp(0.5/abs2(x-S[4]))*exp(0.1/abs2(x-S[5]))*exp(0.5/abs2(x-S[6]))*exp(-0.2x)/sqrt(x)
    singends = [big(0.2);-0.5]
    Itrue = "-3.451882594217500854747246631701525068156038028658631544293802399002537897120929315856041626221574912e-01"
end

# getting parameters
ψ, invψ, dψ, Sψ = Getpsi()
T, β_factor = GetTβ(singends)
dlt, ept = GetStilde(S, Sψ)
m = length(dlt)
D = GetD(dlt)

# calculate C,a,b with NLsolve
print("calculating C,a,b ...\n")
@time C, a, b = Cabsprit(xtoCab(Getparameter().zero))
β = C*β_factor

# nodes n for the trapezoidal rule
if exid == 1
    nidx = [30*i for i = 1:10]
elseif exid == 2 || exid == 3
    nidx = [30*i for i = 1:20]
elseif exid == 4
    nidx = [500*i for i = 1:10]
end

# numerical integration
for i = 1:length(nidx)
    I_1 = nI(HNew, dHNew, nidx[i])
    I_2 = nI(HNew2, dHNew2, nidx[i])
    err_1 = Float16(log10(abs((I_1 - parse(BigFloat,Itrue))/parse(BigFloat,Itrue))))
    err_2 = Float16(log10(abs((I_2 - parse(BigFloat,Itrue))/parse(BigFloat,Itrue))))
    print("n=", nidx[i], " log10|Error|: New1-> ", err_1, "  New2-> ", err_2, "\n")
end
