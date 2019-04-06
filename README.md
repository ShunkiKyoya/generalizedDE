# GeneralizedDE
This is a `Julia 0.6` program for a numerical integration
using a conformal map based on the locations of singularities.
It can be considered to be a generalization of the DE transformations.
We realize high-precision integrations in the case
where the integrand has finite singularities near the real axis.

### Remark
This program uses a package [NLsolve](https://github.com/JuliaNLSolvers/NLsolve.jl)
to solve non-linear systems of equations.

### Settings
You need to assign 4 attributes for integration:
* `f` : the integrand.
* `psid` : the interval.
* `singends` : the degree of endpoint singularities.
* `S` : singularities of the integrand f.

### examples
In `examples.jl`, we show four examples,
where ex.1 and ex.2 are the same as [1].

#### references
[1] R. M. Slevinsky and S. Olver,:
On the use of conformal maps for the acceleration convergence
of the trapezoidal rule and sinc numerical methods, SIAM J. Sci. Comput., 37 (2015),
pp. A676–A700.
