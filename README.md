# GeneralizedDE
This is a `Julia 0.6` program for a numerical integration
using a conformal map based on the locations of singularities.

### Remark
This program uses a package [NLsolve](https://github.com/JuliaNLSolvers/NLsolve.jl)
to solve non-linear systems of equations.

### Settings
You need to assign 4 attributes for integration:
* `f`: the integrand.
* `psid`: the interval and integrand.
* `singends`: the degree of endpoint singularities.
* `S`: singularities of the integrand f.
