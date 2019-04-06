# transformation functions HNew and HNew2

function HNew(t)
    f = C.*sinh.(t-T) + dlt[1]
    for i = 1:m-1
        f = f + 2.*D[i].*atan.(exp.(t-b[i]))
    end
    return f
end

function dHNew(t)
    df = C.*cosh.(t-T)
    for i = 1:m-1
        df = df + D[i]./cosh.(t-b[i])
    end
    return df
end

function HNew2(t)
    f = C.*sinh.(t-T) + (dlt[1]+dlt[m])./2
    for i = 1:m-1
        f = f + D[i].*(big(π)./2.*tanh.((t-b[i]).*2./big(π)))
    end
    return f
end

function dHNew2(t)
    df = C.*cosh.(t-T)
    for i = 1:m-1
        df = df + D[i]./cosh.((t-b[i]).*2./big(π)).^2
    end
    return df
end
