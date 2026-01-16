# Tanh-Sinh Quadrature

The quadrature computes integrals of the form

```math
\mathcal{I}=\int_{-1}^1 f(x) dx.
```

The method is based on a variable transformation which maps the original domain $x \in (-1,1)$ onto the entire real axis $t \in (-\infty, \infty)$ using the transformation:

```math
x = \Psi(t) = \tanh\left(\frac{\pi}{2} \sinh t\right)
```

The derivative (Jacobian) of this transformation is:

```math
\Psi'(t) = \frac{\frac{\pi}{2} \cosh t}{\cosh^2\left(\frac{\pi}{2} \sinh t\right)}
```

The integral then becomes:

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

## Discretization

We approximate the infinite integral using the trapezoidal rule with step size $h$:

```math
\mathcal{I}_h = \sum_{i=-\infty}^{\infty} h g(t_i) = \sum_{i=-\infty}^{\infty} h \Psi'(t_i) f(\Psi(t_i))
```

where $t_i = ih$. Because the transformed integrand $g(t)$ decays double exponentially (like $\exp(-\frac{\pi}{2} e^{|t|})$) as $|t| \to \infty$, we can truncate the infinite sum to a finite window $[-t_n, t_n]$ with negligible error:

```math
\mathcal{I} \approx Q_h^n = \sum_{i=-n}^{n} h \Psi\,'(t_i) \, f(\Psi(t_i))
```

## Error Estimation and Convergence

For an integrand $f(x)$ that is regular in a strip of width $d$ in the complex plane around the interval $[-1, 1]$, the error of the tanh-sinh quadrature decreases exponentially with the number of evaluation points $N = 2n+1$. Specifically, the error is of the order:

```math
|\mathcal{I} - Q_h^n| \approx \mathcal{O}\left(\exp\left(-\frac{\pi d N}{\ln(2 d N)}\right)\right)
```

This rapid convergence rate is the hallmark of double exponential formulas.

## Optimal Step Size and Truncation

The choice of the step size $h$ and the number of points $N$ are coupled. To balance the discretization error (from the trapezoidal rule) and the truncation error (from cutting off the infinite sum), the optimal step size $h$ for a given $N$ is approximately:

```math
h_{opt} \approx \frac{2}{N} \ln(\pi d N)
```

However, in floating-point arithmetic, we are limited by the machine precision. We cannot transform points arbitrarily close to $\pm 1$ without hitting the underflow limit or precision bound of the floating-point type.

This leads to a **maximal step size** constraint to ensure numerical stability. If $t_{max}$ is the largest argument such that we can still distinguish $\Psi(t_{max})$ from $1$ (or weights from $0$), then we must have:

```math
h_{max} = \frac{t_{max}}{n}
```

Typically, $t_{max}$ is determined by the condition where the weights $\Psi'(t)$ underflow to zero or the nodes $\Psi(t)$ become indistinguishable from $\pm 1$ in the given precision.

## Numerical Stability Notes

When dealing with finite precision floating point numbers, numerical instabilities can arise.

The quadrature scheme depends crucially on the evaluations very close to the end-points of the integration domain. The most important cause of numerical instabilities is numerical underflow.

Both the smallest weight and abscissa value are determined by the window size $t_n$. The smallest positive normalized floating point number is $F_{\min}=2^L$.

For the weights to avoid the numerical underflow we need:

```math
t_{\max}^w = \max\{ t\,|\, \Psi\,' (t) \geq F_{\min}\}
```

Similarly the smallest abscissa should exceed the machine epsilon / underflow limit relative to the endpoint:

```math
t_{\max}^x = \max\{ t \,|\, t \leq \Psi^{-1}(1-F_{\min})\}
```

Since these conditions need to be satisfied simultaneously we introduce:

```math
t_{\max}^{xw} = \min\{ t_{\max}^x, \,t_{\max}^w \}
```

It is crucial to ensure $h$ is chosen such that $nh \le t_{max}^{xw}$.
