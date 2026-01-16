# Theory

The **Tanh-Sinh quadrature** (also known as the Double Exponential formula) is a numerical integration method introduced by Hidetosi Takahasi and Masatosi Mori in 1974. It is particularly robust for high-precision integration and functions with end-point singularities.

## Mathematical Formulation

The core idea is to transform a finite integral over $x \in [-1, 1]$ into an infinite integral over $t \in (-\infty, \infty)$ using a variable transformation $x = \tanh(\frac{\pi}{2} \sinh t)$.

Given an integral:

$$I = \int_{-1}^{1} f(x) dx$$

Applying the substitution:

$$x = \tanh\left(\frac{\pi}{2} \sinh t\right)$$

The Jacobian of the transformation (the derivative $dx/dt$) decays double exponentially as $|t| \to \infty$:

$$\frac{dx}{dt} = \frac{\frac{\pi}{2} \cosh t}{\cosh^2\left(\frac{\pi}{2} \sinh t\right)}$$

The integral becomes:

$$I = \int_{-\infty}^{\infty} f\left(\tanh\left(\frac{\pi}{2} \sinh t\right)\right) \cdot \frac{\frac{\pi}{2} \cosh t}{\cosh^2\left(\frac{\pi}{2} \sinh t\right)} dt$$

This transformed integrand decays extremely fast (double exponentially) for large $|t|$. This property allows the infinite sum approximation (trapezoidal rule) to converge very quickly, even if the original function $f(x)$ has singularities at $x = \pm 1$.

## Why Tanh-Sinh?

1.  **Handles Singularities**: The quadrature nodes cluster exponentially densely near $\pm 1$. The weights vanish sufficiently fast to cancel out many algebraic and logarithmic singularities.
2.  **High Precision**: Ideally suited for arbitrary precision arithmetic (e.g., `BigFloat`) where Gaussian quadrature nodes are expensive to compute.
3.  **Simplicity**: It relies on the simple trapezoidal rule on the transformed domain.
