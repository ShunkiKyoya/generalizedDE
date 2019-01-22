# functions for supporting operations

# convert Cab to C,a,b

# convert pole "z" to dlt, ept (by asinh)
function zto_asinh(z)
    dlt0 = real.(asinh.(z))
    ept0 = imag.(asinh.(z))
    id = sortperm(dlt0)
    dlt = [big(dlt0[id[1]])]
    ept = [big(ept0[id[1]])]
    for i = 2:length(dlt0)
        if dlt0[id[i]]-dlt0[id[i-1]]>1e-8
            push!(dlt,dlt0[id[i]])
            push!(ept,ept0[id[i]])
        else
            if ept0[id[i]]<ept[end]
                dlt[end] = dlt0[id[i]]
                ept[end] = ept0[id[i]]
            end
        end
    end
    return dlt,ept
end

# convert pole "z" to dlt, ept (by atanh)
function zto_atanh(z)
    dlt0 = push!(real.(atanh.(z)),big(0.0))
    ept0 = push!(imag.(atanh.(z)),big(π/2))
    id = sortperm(dlt0)
    dlt = [big(dlt0[id[1]])]
    ept = [big(ept0[id[1]])]
    for i = 2:length(dlt0)
        if dlt0[id[i]]-dlt0[id[i-1]]>1e-8
            push!(dlt,dlt0[id[i]])
            push!(ept,ept0[id[i]])
        else
            if ept0[id[i]]<ept[end]
                dlt[end] = dlt0[id[i]]
                ept[end] = ept0[id[i]]
            end
        end
    end
    return dlt,ept
end
# convert pole "z" to dlt, ept (by log)
function zto_log(z)
    dlt0 = real.(log.(z))
    ept0 = imag.(log.(z))
    id = sortperm(dlt0)
    dlt = [big(dlt0[id[1]])]
    ept = [big(ept0[id[1]])]
    for i = 2:length(dlt0)
        if dlt0[id[i]]-dlt0[id[i-1]]>1e-10
            push!(dlt,dlt0[id[i]])
            push!(ept,ept0[id[i]])
        else
            if ept0[id[i]]<ept[end]
                dlt[end] = dlt0[id[i]]
                ept[end] = ept0[id[i]]
            end
        end
    end
    return dlt,ept
end

function zto_logexp(z)
    dlt0 = push!(real.(log.(exp.(z)-1)),big(0.0))
    ept0 = push!(imag.(log.(exp.(z)-1)),big(π))
    id = sortperm(dlt0)
    dlt = [big(dlt0[id[1]])]
    ept = [big(ept0[id[1]])]
    for i = 2:length(dlt0)
        if dlt0[id[i]]-dlt0[id[i-1]]>1e-8
            push!(dlt,dlt0[id[i]])
            push!(ept,ept0[id[i]])
        else
            if ept0[id[i]]<ept[end]
                dlt[end] = dlt0[id[i]]
                ept[end] = ept0[id[i]]
            end
        end
    end
    return dlt,ept
end
