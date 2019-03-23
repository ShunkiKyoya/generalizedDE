# GeneralizedDE
This is a `Julia 0.6` program for a numerical integration
in the case where an integrand has finite singularities.
We construct new transformation formulas which can be seen
as a generalization of DE transformations,
and realize high-precision integrations
in spite of the existence of singularities.

### Remark
This program uses [NLsolve](https://github.com/JuliaNLSolvers/NLsolve.jl)
to solve non-linear systems of equations.
Also, we use [DEQuadrature](https://github.com/MikaelSlevinsky/DEQuadrature.jl)
to compare our method to the previous research[^1].
By default, we use `Plots` and `plotly` to plot.

### Settings
You need to

### References
[^1]:Slevinsky, R. M., Olver, S.,: [On the use of conformal maps for
the acceleration of convergence of the trapezoidal rule and
sinc numerical methods](https://epubs.siam.org/doi/10.1137/140978363),
SIAM J. Sci. Comput., 37(2)(2015),
A676–A700.
