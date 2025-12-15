Approximate the integral such that the error of approximation is less than $10^{-10}$.

1. 
$$\int_0^\infty \frac{1}{1 + 25x^2} dx$$

2.
$$\int_0^1 \frac{\ln(x)}{1 + 25x^2} dx$$

<br></br>

```python
import numpy as np
from scipy import integrate
import mpmath as mp

def report_result(name, result, error, exact=None):
    print(f"{name}")
    print(f"   Numerical result : {result:.15f}")
    print(f"   Estimated error  : {error:.2e}")
    if exact is not None:
        print(f"   Exact value      : {exact:.15f}")
        print(f"   Absolute error   : {abs(result-exact):.2e}")
    print(f"   Requirement met? : {'✓' if error < 1e-10 else '✗'}\n")

def high_precision_check():
    mp.mp.dps = 50

    f = lambda x: mp.log(x) / (1 + 25 * x**2)

    res_mp = mp.quad(f, [0, 1])
    return float(res_mp)

def main():
    print("Numerical integration with required accuracy < 1e-10")
    print("=" * 60)

    f1 = lambda x: 1 / (1 + 25*x**2)
    res1, err1 = integrate.quad(f1, 0, np.inf, epsabs=1e-12, epsrel=1e-12)
    report_result("1. ∫₀∞ 1/(1+25x²) dx", res1, err1, np.pi/10)

    f2 = lambda x: np.log(x)/(1+25*x**2)
    res2, err2 = integrate.quad(f2, 0, 1, epsabs=1e-12, epsrel=1e-12, points=[0])
    report_result("2. ∫₀¹ ln(x)/(1+25x²) dx", res2, err2)

    res2_mp = high_precision_check()
    print("High-precision validation with mpmath (50 digits):")
    print(f"   Result (mpmath)  : {res2_mp:.15f}")
    print(f"   Absolute error   : {abs(res2 - res2_mp):.2e}")
    print()

    print("=" * 60)
    print("Final Results:")
    print(f"∫₀∞ 1/(1+25x²) dx      = {res1:.15f}")
    print(f"∫₀¹ ln(x)/(1+25x²) dx   = {res2:.15f}")

if __name__ == "__main__":
    main()

    print("\nSummary:")
    print("1. Used scipy.integrate.quad to compute both integrals with error < 1e-10.")
    print("2. Verified the first integral against its exact value π/10.")
    print("3. Verified the second integral using mpmath high-precision integration.")
```
#### Conclusion (Numerical Integration with Required Accuracy $$< 10^{-10}$$)
 

**Integration Results**

* $$\int_0^\infty \frac{1}{1+25x^2} dx$$

| Metric | Value |
| :---: | :---:|
| Numerical Result | $0.314159265358979$ |
| Estimated Error | $2.72 \times 10^{-13}$ |
| Exact Value | $0.314159265358979$ |
| Absolute Error | $0.00 \times 10^{+00}$ |
| Requirement Met? | yes |

<br></br>
* $$\int_0^1 \frac{\ln(x)}{1+25x^2} dx$$

| Metric | Value |
| :---: | :---: |
| Numerical Result | $$-0.545444563419798$$ |
| Estimated Error | $$7.65 \times 10^{-13}$$ |
| Requirement Met? | yes |

<br></br>
* High-precision validation with mpmath (50 digits):

| Metric | Value |
| :---: | :---: |
| **Result (mpmath)** | $-0.545444563419798$ |
| **Absolute Error** | $2.22 \times 10^{-16}$ |

<br></br>

* Final Summary

| Integral | Result |
| :---: | :---: |
| $\int_0^\infty \frac{1}{1+25x^2} dx$ | $0.314159265358979$ |
| $\int_0^1 \frac{\ln(x)}{1+25x^2} dx$ | $-0.545444563419798$ |

  * Used `scipy.integrate.quad` to compute both integrals, successfully achieving an estimated error less than the required tolerance of $10^{-10}$.
  * The first integral, $\int_0^\infty \frac{1}{1+25x^2} dx$, was analytically verified against its exact value, which is $\frac{\pi}{10}$.
  * The second integral, $\int_0^1 \frac{\ln(x)}{1+25x^2} dx$, was verified using the high-precision integration capabilities of the `mpmath` library.

<br></br>

**Programing**
link:https://colab.research.google.com/drive/172BRKa3BMVr5Pu2PM1V2m4YW264hcnAc?usp=sharing

