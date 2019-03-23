function Getpsi()
    if psid == 1
        ψ = x -> tanh.(x)
        invψ = x -> atanh.(x)
        dψ = x -> 1./cosh.(x).^2
        Sψ = [complex(big(0.0),big(π)/2)]
    elseif psid == 2
        ψ = x -> sinh.(x)
        invψ = x -> asinh.(x)
        dψ = x -> cosh.(x)
        Sψ = []
    elseif psid == 3
        ψ = x -> exp.(x)
        invψ = x -> log.(x)
        dψ = x -> exp.(x)
        Sψ = []
    elseif psid == 4
        ψ = x -> log.(exp.(x)+1)
        invψ = x -> log.(exp.(x)-1)
        dψ = x -> exp.(x)./(exp.(x)+1)
        Sψ = [complex(big(0.0),big(π))]
    end
    return ψ, invψ, dψ, Sψ
end

function GetTβ(singends)
    if psid == 1
        T = log((singends[1]+1)/(singends[2]+1))/2
        β_factor = sqrt((singends[1]+1)*(singends[2]+1))
    elseif psid == 2
        T = log((singends[1]+1)/(singends[2]+1))/2
        β_factor = sqrt((singends[1]+1)*(singends[2]+1))/2
    elseif psid == 3
        T = log(-(singends[1]+1)/(singends[2]+1))/2
        β_factor = sqrt(-(singends[1]+1)*(singends[2]+1))/2
    elseif psid == 4
        T = log(singends[1]/(singends[2]+1))/2
        β_factor = sqrt(singends[1]*(singends[2]+1))/2
    end
    return T, β_factor
end

function GetStilde(S,Sψ)
    if length(Sψ) == 1
        dlt0 = push!(real.(invψ(S)),real(Sψ[1]))
        ept0 = push!(imag.(invψ(S)),imag(Sψ[1]))
    else
        dlt0 = real.(invψ(S))
        ept0 = imag.(invψ(S))
    end
    id = sortperm(dlt0)
    dlt = [big(dlt0[id[1]])]
    ept = [big(ept0[id[1]])]
    for i = 2:length(dlt0)
        if dlt0[id[i]] != dlt0[id[i-1]]
            push!(dlt,dlt0[id[i]])
            push!(ept,ept0[id[i]])
        elseif ept0[id[i]]<ept[end]
                dlt[end] = dlt0[id[i]]
                ept[end] = ept0[id[i]]
        end
    end
    return dlt,ept
end

function GetD(dlt)
    D = zeros(BigFloat,m-1)
    for i = 1:m-1
        D[i] = (dlt[i+1]-dlt[i])/big(π)
    end
    return D
end
