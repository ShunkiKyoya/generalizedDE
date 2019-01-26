# "psid" defines psi.
# psid = 1: ψ = tanh(.)
#        2: ψ = sinh(.)
#        3: ψ = exp(.)
#        4: ψ = log(exp(.)+1)
#
# "singends" defines degree of endpoint singularities.
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
#
#
#

using Base.Test
import Base.digits

### This program uses NLsolve for solving nonliniar equations.
using NLsolve

### setting precision of BigFloat
digits(n::Int) = setprecision(round(Int,ceil(n*log2(10))))
digits(300) # significant digits

include("Getparameter.jl")
include("Trapezoidal.jl")
include("psiandT.jl")
include("HNew.jl")

exid = 2
#なぜか例2が動かない

if exid == 1
    psid = 1
    S = [complex(big(-0.5),1.0); complex(0.5,0.5)]
    f = x-> exp.(1./abs2.(x-S[1])).*log.(1-x)./abs2.(x-S[2])./sqrt.(1+x)
    Itrue = "-2.04645081160694748690442050179886173463698400851312978159495108281833923937599915411239665314732909414150138401559892277212580897054853305044117124435121918549124604645789190248721167573993297627891884242759929922470299331035480062586543433278436544374899377193299256405496172603459957078880800178958528771642962893695508317394167405216487683527096100475182610564645687173201068270992551244010304597e+00"
    singends = [big(0.0); -0.5]
elseif exid == 2
    psid = 2
    S = [complex(big(-2.0),1.0);complex(-1.0,0.5);complex(1.0,0.25);complex(2.0,1.0)]
    f = x-> exp.(10./abs2.(x-S[1])).*cos.(10./abs2.(x-S[2]))./abs2.(x-S[3])./abs.(x-S[4])
    Itrue = "1.50133619876062770101030470326173553208854739646240081225845195322624377330867094144192464692844931351956249439446687259575908821074467954700173593866783934590906808530844595015500541680530415264580194453133658724574536164866695611699188173914168980057456962540582381978462005370937587010012338985572909021270993977527355156445163416743376725727256310809095229447802636632237845095430431700882602931e+01"
    singends = [-big(3.0);-3.0]
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

INew = zeros(BigFloat,2)
nidx = [100;200]
nIarray!(INew,HNew,dHNew,nidx)
