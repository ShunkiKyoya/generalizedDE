# functions for numerical integration with trapezoidal rule

function nI(H,dH,n)
    d = big(π)/2
    F = t -> f(ψ(H(t))).*dψ(H(t)).*dH(t)
    I = big(0.0)
    h = log(big(π)*d*n*2/β)/n
    for j = -n:n
        if !isnan(F(j*h)) && !isinf(F(j*h)) && !isinf(-F(j*h))
            I = I + F(j*h)
        end
    end
    I = I*h
    return I
end
