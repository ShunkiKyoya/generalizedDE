# functions for getting parameters of the proposed transformation formula

function Cabsprit(Cab)
    C = Cab[1]
    a = zeros(BigFloat,m)
    b = zeros(BigFloat,m-1)
    for i = 1:m-1
        a[i] = Cab[2*i]
        b[i] = Cab[2*i+1]
    end
    a[m] = Cab[2*m]
    return C,a,b
end

# convert x to Cab
function xtoCab(x)
    Cab = zeros(BigFloat,2*m)
    Cab[1] = exp(x[1])
    for i = 2:2*m
        Cab[i] = Cab[i] + x[2]
    end
    for i = 3:2*m
        for j = i:2*m
            Cab[j] = Cab[j] + exp(x[i])
        end
    end
    return Cab
end

# jacobian of xtoCab
function xtoCab_Jacobi(x)
    J = zeros(BigFloat,2*m,2*m)
    J[1,1] = exp(x[1])
    for i = 2:2*m
        J[i,2] = 1
    end
    for j = 3:2*m
        for i = j:2*m
            J[i,j] = exp(x[j])
        end
    end
    return J
end

# 2*Np equations
function mEq_Cab(Cab)
    F = zeros(BigFloat,2*m)
    C = Cabsprit(Cab)[1]
    a = Cabsprit(Cab)[2]
    b = Cabsprit(Cab)[3]

    for i = 1:m
        F[i] = C*sinh(a[i]-T)
        F[i+m] = C*cosh(a[i]-T) - ept[i]
        for j = 1:m-1
            F[i] = F[i] - D[j]/sinh(a[i]-b[j])
            if i>j
                F[i+m] = F[i+m] - D[j]*log(tanh((a[i]-b[j])/2))
            else
                F[i+m] = F[i+m] - D[j]*log(tanh((b[j]-a[i])/2))
            end
        end
    end
    return F
end

# jacobian of "NpEq_Cab"
function mEq_Cab_Jacobi(Cab)
    J = zeros(BigFloat,2*m,2*m)
    C = Cabsprit(Cab)[1]
    a = Cabsprit(Cab)[2]
    b = Cabsprit(Cab)[3]

    for i = 1:m
        J[i   ,1] = sinh(a[i]-T)
        J[i+m,1] = cosh(a[i]-T)
        J[i   ,2*i] = C*cosh(a[i]-T)
        J[i+m,2*i] = C*sinh(a[i]-T)
        for j = 1:m-1
            J[i   ,2*i] = J[i,2*i] + D[j]*cosh(a[i]-b[j])/sinh(a[i]-b[j])^2
            J[i+m,2*i] = J[m+i,2*i] - D[j]/sinh(a[i]-b[j])
            J[i   ,2*j+1] = - D[j]*cosh(a[i]-b[j])/sinh(a[i]-b[j])^2
            J[i+m,2*j+1] = D[j]/sinh(a[i]-b[j])
        end
    end
    return J
end

function mEq_x!(F,x)
    F .= mEq_Cab(xtoCab(x))
end

function mEq_x_Jacobi!(J,x)
    J .= mEq_Cab_Jacobi(xtoCab(x))*xtoCab_Jacobi(x)
end

function Getparameter()
    myans = nlsolve(mEq_x!,mEq_x_Jacobi!,zeros(BigFloat,2*m))
    @test converged(myans)
    return myans
end
