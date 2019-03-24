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
Also, it uses [DEQuadrature](https://github.com/MikaelSlevinsky/DEQuadrature.jl)
to compare our method to the previous research [1].
By default, it uses `Plots` and `plotly` to plot.

### Settings
* `psid` defines the interval and integrand:

| psid | interval | integrand |
|:---:|:---:|:---:|
|1|(-1,1)|f(x)|
|2|(-∞, ∞)|f(x)|
|3|(0,∞)|f(x)|
|4|(0,∞)|f1(x)exp(-vx)|
### References
[1] Slevinsky, R. M., Olver, S. (2005) :
[On the use of conformal maps for the acceleration of convergence of the trapezoidal rule
and sinc numerical methods](https://epubs.siam.org/doi/10.1137/140978363),
SIAM J. Sci. Comput., 37(2), A676–A700.
