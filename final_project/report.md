**About**

This paper addresses the problem of interpolating real or complex data sampled at equispaced points on the interval $$[-1,1]$$.

Core idea of the AAA method

The AAA algorithm computes a rational approximation to the given data with a default relative tolerance of $$10^{-13}$$. In MATLAB (Chebfun), the basic command is


$$r = aaa(F)$$


If a polynomial rather than a rational function is required, the rational approximant can be converted afterward.

**Numerical interpolant vs. exact interpolation**

The output of AAA is generally not a mathematically exact interpolant (i.e., it does not necessarily have degree $$n-1$$). Instead, it produces a numerical interpolant. This distinction is crucial for robustness.

Example (Fig. 1)

Consider

$$
f(x) = \frac{e^x}{\sqrt{1+9x^2}}
$$

sampled at 50 equispaced points. AAA produces a degree-17 rational function with an error of only

$$
9.6 \times 10^{-14}.
$$

By contrast, the exact degree-49 polynomial interpolant suffers from severe Runge phenomena, with an error as large as

$$
109.3.
$$

**Accuracy–stability tradeoff**

For equispaced data, there is a fundamental tradeoff between accuracy and stability. An impossibility theorem states that

* Exponential convergence typically comes with exponential instability.
* Stable algorithms can achieve at best root-exponential convergence, $$exp(-C\sqrt{n}).$$


<br></br>

**Existing Methods**

The paper compares AAA against five major competing approaches.

* Polynomial Least Squares (PLS)

Direct polynomial interpolation with degree $$d=n-1$$ leads to error growth of order

$$
O(2^n).
$$

Using a reduced degree $$d \approx n/\gamma$$ with oversampling factor $$\gamma>1$$ mitigates this instability. For example, when $$\gamma=2$$, the growth rate is reduced to approximately

$$
(1.14)^n.
$$

* Fourier Extension

The function is approximated by a Fourier series defined on an extended interval $$[-T,T]$$. In the experiments, $$T=2$$ , oversampling factor $$=2.$$ are used.

* Fourier Series with Corrections

This method combines a trigonometric approximation with polynomial correction terms to handle endpoint discontinuities. In the experiments, the Fourier series is augmented by a polynomial of degree approximately $$\sqrt{n}$$.

* Splines

Cubic splines are employed. They are extremely stable but converge only algebraically, with a fixed rate of

$$
O(n^{-4}).
$$

* Floater–Hormann Rational Interpolation (Chebfun `equi`)

This is a barycentric rational interpolant guaranteed to have no poles in $$[-1,1]$$. Chebfun’s `equi` option adaptively selects the interpolation degree.

<br></br>

**Numerical Comparison**



AAA is compared with the five methods above using the following five test functions:

* $$f_A(x) = \sqrt{1.21 - x^2}$$ (branch points at $$\pm1.1$$)
* $$f_B(x) = \sqrt{0.01 + x^2}$$ (branch points at $$\pm0.1i$$, very close to the real axis)
* $$f_C(x) = \tanh(5x)$$ (poles on the imaginary axis)
* $$f_D(x) = \sin(40x)$$ (entire function, highly oscillatory)
* $$f_E(x) = \exp(-1/x^2)$$ ($$C^\infty$$ but not analytic)

Main results (Fig. 2)

* Dominance of AAA:AAA outperforms all other methods in almost every case and is consistently the first to reach an accuracy of $$10^{-10}$$.

* Capturing singularities: For $$f_C(x)=\tanh(5x)$$, AAA rapidly detects the nearby poles, leading to extremely fast convergence.

* Numerical interpolation behavior:For small $$n$$, AAA behaves like an exact interpolant. As $$n$$ increases, it transitions into a numerical interpolant: the error at grid points remains around $$10^{-13}$$ while avoiding overfitting.

* Handling spurious poles (AAA-LS):If AAA produces undesirable poles inside $$[-1,1]$$, the paper proposes AAA–Least Squares (AAA-LS):

1. Remove the bad poles.
2. Retain the remaining poles as a basis.
3. Recompute the approximation via least squares.


<br></br>
**Convergence Properties**

**Typical AAA convergence phases**

Initial phase (small $$n$$) : Resolution is insufficient; spurious poles may appear, triggering a switch to least squares.Rapid convergence phase (middle phase): This is where AAA excels. It exploits the analytic structure (poles, branch points). Grid-point errors are around $$10^{-13}$$, while inter-grid errors may still be larger (e.g., $$10^{-6}$$). As $$n$$ increases, the spacing shrinks and the global error rapidly decreases.Leveling-off phase: The error plateaus due to the prescribed tolerance of $$10^{-13}$$.


The Amber function test (Figs. 3–4)

The authors construct an "Amber function" $$A(x)$$ whose coefficients are determined by the binary digits of $$\pi$$. It is analytic inside a Bernstein 2-ellipse but has essentially no exploitable structure.

Result: For this nearly non-rational function, AAA does not outperform the Floater–Hormann method.

AAA’s advantage lies in its ability to exploit analytic structure. AAA is a nonlinear method. Unlike linear schemes that always use all $$n$$ degrees of freedom, AAA automatically stops increasing the rational degree once the target accuracy is reached—effectively choosing an optimal oversampling rate. The paper reveals an unexpected phenomenon: the stability of Fourier extensions arises from the use of an ill-conditioned basis (complex exponentials $$e^{i\pi kx/2}$$) combined with floating-point roundoff. These roundoff errors limit the growth of the condition number and thus accidentally stabilize the computation.
If a well-conditioned basis (e.g., Arnoldi) is used instead, Fourier extensions become exponentially unstable.

<br></br>

**Discussion**

Analytic continuation: AAA not only provides high accuracy on $$[-1,1]$$ but also yields excellent analytic continuation into the complex plane (Fig. 8).

Missing data: AAA can handle equispaced grids with missing values.

Computational complexity: The computational cost is $$O(md^3),$$ where $$m$$ is the number of data points and $$d$$ is the rational degree. While very large degrees ($$d>100$$) can be expensive, in practice $$d$$ is usually small.

Noisy data: For noisy data, the tolerance `tol` should be set one or two orders of magnitude above the noise level.

<br></br>


**Summary**

The core mathematical message of the paper lies in comparing error growth rates—such as $$2^n$$ versus $$C^n$$—and in demonstrating how AAA, through a barycentric rational representation

$$
r(x) = \frac{p(x)}{q(x)},
$$

combined with a nonlinear, adaptive algorithm, achieves a breakthrough in the classical problem of interpolation on equispaced grids.
