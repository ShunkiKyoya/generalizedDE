# functions for getting parameters of the proposed transformation formula

function Cabsprit(Cab)
    C = Cab[1]
    a = zeros(BigFloat,Np)
    b = zeros(BigFloat,Np-1)
    for i = 1:Np-1
        a[i] = Cab[2*i]
        b[i] = Cab[2*i+1]
    end
    a[Np] = Cab[2*Np]
    return C,a,b
end

# convert x to Cab
function xtoCab(x)
    Cab = zeros(BigFloat,2*Np)
    Cab[1] = exp(x[1])
    for i = 2:2*Np
        Cab[i] = Cab[i] + x[2]
    end
    for i = 3:2*Np
        for j = i:2*Np
            Cab[j] = Cab[j] + exp(x[i])
        end
    end
    return Cab
end

# jacobian of xtoCab
function xtoCab_Jacobi(x)
    J = zeros(BigFloat,2*Np,2*Np)
    J[1,1] = exp(x[1])
    for i = 2:2*Np
        J[i,2] = 1
    end
    for j = 3:2*Np
        for i = j:2*Np
            J[i,j] = exp(x[j])
        end
    end
    return J
end

# 2*Np equations
function NpEq_Cab(Cab)
    F = zeros(BigFloat,2*Np)
    C = Cabsprit(Cab)[1]
    a = Cabsprit(Cab)[2]
    b = Cabsprit(Cab)[3]

    for i = 1:Np
        F[i] = C*sinh(a[i]-T)
        F[i+Np] = C*cosh(a[i]-T) - ept[i]
        for j = 1:Np-1
            F[i] = F[i] - D[j]/sinh(a[i]-b[j])
            if i>j
                F[i+Np] = F[i+Np] - D[j]*log(tanh((a[i]-b[j])/2))
            else
                F[i+Np] = F[i+Np] - D[j]*log(tanh((b[j]-a[i])/2))
            end
        end
    end
    return F
end

# jacobian of "NpEq_Cab"
function NpEq_Cab_Jacobi(Cab)
    J = zeros(BigFloat,2*Np,2*Np)
    C = Cabsprit(Cab)[1]
    a = Cabsprit(Cab)[2]
    b = Cabsprit(Cab)[3]

    for i = 1:Np
        J[i   ,1] = sinh(a[i]-T)
        J[i+Np,1] = cosh(a[i]-T)
        J[i   ,2*i] = C*cosh(a[i]-T)
        J[i+Np,2*i] = C*sinh(a[i]-T)
        for j = 1:Np-1
            J[i   ,2*i] = J[i,2*i] + D[j]*cosh(a[i]-b[j])/sinh(a[i]-b[j])^2
            J[i+Np,2*i] = J[Np+i,2*i] - D[j]/sinh(a[i]-b[j])
            J[i   ,2*j+1] = - D[j]*cosh(a[i]-b[j])/sinh(a[i]-b[j])^2
            J[i+Np,2*j+1] = D[j]/sinh(a[i]-b[j])
        end
    end
    return J
end

function NpEq_x!(F,x)
    F .= NpEq_Cab(xtoCab(x))
end

function NpEq_x_Jacobi!(J,x)
    J .= NpEq_Cab_Jacobi(xtoCab(x))*xtoCab_Jacobi(x)
end

function Getparameter()
    myans = nlsolve(NpEq_x!,NpEq_x_Jacobi!,zeros(BigFloat,2*Np))
    @test converged(myans)
    return myans
end
