# functions for numerical integration with trapezoidal rule

function nI(H,dH,n)
    d = big(π)/2
    F = t -> f(ψ(H(t))).*dψ(H(t)).*dH(t)
    I = big(0.0)
    h = log(big(π)*d*n*2/β)/n
    print("n = ", n, " (trapezoidal rule) \n")
    @time for j = -n:n
        if !isnan(F(j*h)) && !isinf(F(j*h)) && !isinf(-F(j*h))
            I = I + F(j*h)
        end
    end
    I = I*h
    return I
end

function nIarray!(Iarray,H,dH,nidx)
    for i = 1:length(Iarray)
        Iarray[i] = nI(H,dH,nidx[i])
    end
end
