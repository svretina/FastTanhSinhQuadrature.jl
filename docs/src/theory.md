# Tanh-Sinh Quadrature

The quadrature computes integrals of the form

```math
\mathcal{I}=\int_{-1}^1 f(x) dx.
```

The method is based on a variable transformation which maps the original domain $x \in (-1,1)$ onto the entire real axis $t \in (-\infty, \infty)$.

```math
\mathcal{I} = \int_{-\infty}^\infty g(t) dt, \quad g(t)= f(\Psi(t)) \, \Psi^\prime (t)
```

Since the method is specific for the $x \in (-1,1)$ domain, one must cast the desired integral to this domain by the linear substitution:

```math
x(u)=\frac{b+a}{2}+\frac{b-a}{2}u
```

This transformation changes an arbitrary interval $[a,b]$ to $[-1,1]$, hence

```math
\int_a^b f(x)dx= \frac{b-a}{2}\int_{-1}^1 f(x(u))du = \frac{b-a}{2}\int_{-\infty}^\infty f(x(u(t))) w(t) dt
```

![Figure 1: Transformation visualization (Source: arXiv:2007.15057)](figure_1.pdf)

Now we need to integrate the whole real axis numerically. The Trapezoidal rule is the most efficient among quadratures with equidistant abscissae. Applying the trapezoidal rule for a step size $h$:

```math
\mathcal{I}_h= \lim_{m\rightarrow \infty} \sum_{i=-m}^{m} h \Psi\,'(t_i) \, f(\Psi(t_i))
```

Keeping only $2n+1$ function evaluations:

```math
\mathcal{I}_h= \sum_{i=-n}^{n} h \Psi\,'(t_i) \, f(\Psi(t_i))\quad \text{with} \quad t_i \in [-t_n,t_n]
```

In the tanh-sinh quadrature the window size $[-t_n,t_n]$ is not fixed *a priori*.

## Numerical Stability

When dealing with finite precision floating point numbers, numerical instabilities can arise.

The quadrature scheme depends crucially on the evaluations very close to the end-points of the integration domain, where the abscissae are very high (close to $\pm 1$).

The most important cause of numerical instabilities is numerical underflow, occurring when numerical values become smaller than the underflow level.

Both the smallest weight and abscissa value are determined by the window size $t_n$. The smallest positive normalized floating point number is $F_{\min}=2^L$ where $L$ is the smallest exponent representable in a given floating point model (e.g., `floatmin(T)` in Julia).

For the weights to avoid the numerical underflow we need:

```math
t_{\max}^w = \max\{ t\,|\, \Psi\,' (t) \geq F_{\min}\}
```

In code: `w ≥ floatmin(T)`

Similarly the smallest abscissa should exceed the machine epsilon / underflow limit relative to the endpoint:

```math
t_{\max}^x = \max\{ t \,|\, t \leq \Psi^{-1}(1-F_{\min})\}
```

Since these conditions need to be satisfied simultaneously we introduce:

```math
t_{\max}^{xw} = \min\{ t_{\max}^x, \,t_{\max}^w \}
```

## Notes

The condition $t \ge \Psi^{-1}(1-\rm{F_{\min}})$ is to avoid $f(\Psi(t)) = \pm 1$ where singularities can lie. `eps(T)` should generally be used instead of `floatmin(T)` for checking distance from 1.

When transforming a general domain $[a,b]$ to $[-1,1]$ we introduce further floating point operations that are not accounted into the original conditions to avoid underflow. One strictly must apply the condition:

```math
t \ge \Psi^{-1}((1-\rm{F_{\min}}-(a+b)/2)/(b-a)/2)
```

This though defeats the whole purpose of our design of a reusable object since one would need to pass the domain of integration at each step and create the points from scratch.

This problem can be solved if we go to a slightly different transformation that uses higher odd powers of $t$ in the main transformation (Double Exponential formulas with higher order decay).

## Bibliography

1. [Tanh-sinh quadrature for single and multiple integration using floating-point arithmetic](https://arxiv.org/abs/2007.15057)
2. [Efficient computation of Sommerfeld integral tails-methods and algorithms](https://www.researchgate.net/publication/301544502_Efficient_computation_of_Sommerfeld_integral_tails-methods_and_algorithms)
3. [DoubleExponentialFormulas.jl](https://machakann.github.io/DoubleExponentialFormulas.jl/stable/)
