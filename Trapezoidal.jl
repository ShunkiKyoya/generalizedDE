# functions for numerical integration with trapezoidal rule

function nI(F,n,β,d)
    I = big(0.0)
    l = log(big(π)*d*n*2/β)/n
    print(n, "time\n")
    @time for j = -n:n
        if !isnan(F(j*l)) && !isinf(F(j*l)) && !isinf(-F(j*l))
            I = I + F(j*l)
        end
    end
    I = I*l
    return I
end

function nIarray!(Iarray,F,nidx,β,d)
    for i = 1:length(Iarray)
        Iarray[i] = nI(F,nidx[i],β,d)
    end
end
