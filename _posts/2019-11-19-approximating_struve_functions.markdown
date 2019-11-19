---
layout: post
title: Approximating Struve functions
category: posts
draft: false
---

As special functions go, it seems the Struve function goes somewhat unloved. Whilst one cannot move for implementations of Bessel function, both the default MATLAB package and the Julia programming language are left lacking, leaving the programmer to roll their own.

Sadly, the usual trick of heading over to Wikipedia to look an easy approximation doesn't end with easy answers. Likewise, whilst Zhang and Jin has a chapter on Struve functions, it doesn't go beyond the power series expansion found in Abramowitz and Stegan. 

This left me a little lost. Olver et al. suggest the direct evaluation of power series, but don't give guidance on the number of terms or when to switch to the large z expansion. Other authors present their approximations, but don't delve extensively into the error of their approximation. Here, I aim to explore this a little.

# A little about Struve functions

Struve functions are a component of a general solution to a nasty differential equation. However, for me they showed up in integral form. As a special function, there is no neat means by which to evaluate them. The Struve functions are $$\mathbf{H}_v(z)$$ and $$\mathbf{K}_v(z)$$, and the modified Struve functions are $$\mathbf{L}_v(z)$$ and $$\mathbf{M}_v(z)$$. $$\mathbf{K}_v(z)$$ and $$\mathbf{M}_v(z)$$ can be obtained from $$\mathbf{H}_v(z)$$ and $$\mathbf{L}_v(z)$$ through the subtraction of Bessel functions and modified Bessel functions. Since approximations of Bessel functions are everywhere, I consider this a non-issue. Likewise, the Struve function $$\mathbf{H}_v(z)$$ is related to the modified Struve function $$\mathbf{L}_v(z)$$ by the formula

$$\mathbf{L}_v(z) = - i e^{-\frac{iv\pi}{2}} \mathbf{H}_v(iz)$$

meaning that so long as Struve function $$\mathbf{H}_v(iz)$$ can be evaluated, we can obtain every other variety.

Furthermore, recurrence relations can be used to obtain Struve functions of any integer order given $$\mathbf{H}_0(z)$$ and $$\mathbf{H}_1(z)$$.

$$\mathbf{H}_{v-1}(z) + \mathbf{H}_{v+1}(z) = \frac{2v}{z}\mathbf{H}_{v}(z) + \frac{\left(\frac{1}{2}z\right)^v}{\sqrt{\pi}\Gamma\left(v + \frac{3}{2}\right)}$$

Consequently, this post focuses on $$\mathbf{H}_0(z)$$ and $$\mathbf{H}_1(z)$$, seemingly the most studied orders. I won't consider non-integer orders.

# Power series expansion

Looking in Abramowitz and Stegan (or Olver et al.) will give you a definition of the Struve function $$\mathbf{H}_v(z)$$. However, I'll jump straight to the most obvious method of evaluation - the power series expansion:

$$
\mathbf{H}_v(z) = \left(\frac{1}{2}z\right)^{v+1}
\sum_{k=0}^{\infty} \frac{(-1)^k\left(\frac{1}{2}z\right)^{2k}}{\Gamma(k+\frac{3}{2}) \Gamma(k + v + \frac{3}{2})}
$$

and for $$v=0$$ and $$v=1$$

$$
\mathbf{H}_0(z) = \frac{2}{\pi}
\left( z - \frac{z^3}{1^2\cdot3^2} + \frac{z^5}{1^2\cdot3^2\cdot5^2} - ... \right)  
$$

$$
\mathbf{H}_1(z) = \frac{2}{\pi}
\left( \frac{z^2}{1^2\cdot3} - \frac{z^4}{1^2\cdot3^2\cdot5}  + \frac{z^6}{1^2\cdot3^2\cdot5^2\cdot7} -... \right)  
$$

These series are easy to program and always convergent. Care must be taken to using floating point in the denominator - an integer will rapidly overflow. 

A large number of terms are required to achieve good accuracy. To obtain a solution with a relative error of approximately 1e-6, the number of terms required in the $$\mathbf{H}_0(z)$$ expansion are plotted:

![H_0(z) number of terms for 1e-6 convergence]({{ site.url }}/images/2019-11-19-approximating-struve-functions/num_terms_H0_polynomial.png "H_0(z) number of terms for 1e-6 convergence")

And for $$\mathbf{H}_1(z)$$:

![H_1(z) number of terms for 1e-6 convergence]({{ site.url }}/images/2019-11-19-approximating-struve-functions/num_terms_H1_polynomial.png "H_1(z) number of terms for 1e-6 convergence")

The number of terms required gets bigger scales with $$\vert z \vert$$.

For large values of $$z$$ Abramowitz and Stegan give asymptotic expansions:

$$
\mathbf{H}_0(z) = \frac{2}{\pi}
\left( \frac{1}{z} - \frac{1}{z^3} + \frac{1^2\cdot3^2}{z^5} - \frac{1^2\cdot3^2\cdot5^2}{z^7} + ... \right) + Y_0(z)
$$

$$
\mathbf{H}_1(z) = \frac{2}{\pi}
\left( 1 + \frac{1}{z^2} - \frac{1^2\cdot3}{z^4} + \frac{1^2\cdot3^2\cdot5}{z^6} + ... \right) + Y_1(z)
$$

The numerator in each term of these expansions grows faster than the denominator, to my mind implying infinite truncation error for any finite argument. None the less, taking the just the first term of these and combining with the large value expansion $$Y_v(z) = \sqrt{2 / (\pi z)}\sin(z - \pi v / 2- \pi / 4)$$ (A&S 9.2.2) gives:

$$
\mathbf{H}_0(z) \approx \frac{2}{\pi z} + \sqrt{\frac{2}{\pi z}} \sin\left(z - \frac{\pi}{4}\right) \text{ as }\vert z\vert \rightarrow \infty
$$

$$
\mathbf{H}_1(z) \approx \frac{2}{\pi} + \sqrt{\frac{2}{\pi z}} \sin\left(z - \frac{3\pi}{4}\right) \text{ as }\vert z\vert \rightarrow \infty
$$

Naturally, its unclear what constitutes large $$z$$. I therefore compared these approximations to the Python mpmath implementation using 25 digits of precisions for real x. The relative accuracy of these approximations was disappointing:

![Struve large x convergence 1 term]({{ site.url }}/images/2019-11-19-approximating-struve-functions/big_x_approximation_1_term.png "Struve large x convergence 1 term")

The main source of error here appears to be the Bessel function approximation. Using instead

$$
\mathbf{H}_0(z) \approx \frac{2}{\pi z} + Y_0(z)
$$

$$
\mathbf{H}_1(z) \approx \frac{2}{\pi} + Y_1(z)
$$

Additionally we can add additional terms to the expansion. The results for real x  are show:

![Struve large x convergence using better Bessel]({{ site.url }}/images/2019-11-19-approximating-struve-functions/big_x_approximation_exact_bessel_multiterm.png "Struve large x convergence using better Bessel")

Certainly the approximation of the Bessel function results in significant errors. Presumably the noise in the $$\mathbf{H}_0(x)$$ 1 term curve and the 5 term curves is due to the introduction of double precision numbers into the calculations. The addition of multiple terms is also very beneficial. The implications of computing the approximation in double precision is noticeable.

A composite of the two methods can be used. If a 60 term power series is used for $$\vert z \vert < 20$$ and the large value expansion with 7 terms is used otherwise, the following results are obtained for $$\mathbf{H}_0(z)$$:

![Composite power series accuracy H0]({{ site.url }}/images/2019-11-19-approximating-struve-functions/composite_h0.png "Composite power series accuracy H0")

And for $$\mathbf{H}_1(z)$$:
![Composite power series accuracy H1]({{ site.url }}/images/2019-11-19-approximating-struve-functions/composite_h1.png "Composite power series accuracy H1")

Reasonable results can be obtained. However, for large portions of the domain, the computational cost could be reduced significantly by changing the number of terms. The most difficult region seems to be when the real part is around 20. Accuracy is only single precision despite the computational effort. 

# Aarts and Janssen's approximations

Aarts & Janssen develop an approximation as follows:

$$
\mathbf{H}_1(z) \approx \frac{2}{\pi} - J_0(z) + \left(\frac{16}{\pi} - 5\right)\frac{\sin(z)}{z} + \left(12 - \frac{36}{\pi}\right)\frac{1 - \cos(z)}{z^2}
$$

The log10 of the relative error can be plotted.

![Aarts & Janssen h1 error]({{ site.url }}/images/2019-11-19-approximating-struve-functions/aarts_h1_err.png "Aarts & Janssen h1 error")

and $$\mathbf{H}_0(z)$$ as

$$
\mathbf{H}_0(z) \approx J_1(z) + \left(7-\frac{20}{\pi}\right)\frac{1-\cos(z)}{z} + \left(\frac{36}{\pi} - 12\right)\frac{\sin(z) - z\cos(z)}{z^2}
$$

![Aarts & Janssen h0 error]({{ site.url }}/images/2019-11-19-approximating-struve-functions/aarts_h0_err.png "Aarts & Janssen h1 error")

Evidently, whilst these expressions are convenient, they leave something to be desired with respect to accuracy.


# Piecewise polynomial approximation

Newman describes a composite approach for real $$x$$.

For real $$x$$ where $$0 \leq x \leq 3$$:

$$
\mathbf{H}_0(x) =  1.909859164 \left(\frac{x}{3}\right) - 1.909855001 \left(\frac{x}{3}\right)^3\\+ 0.687514637 \left(\frac{x}{3}\right)^5-0.126164557 \left(\frac{x}{3}\right)^7\\+ 0.013828813 \left(\frac{x}{3}\right)^9-0.000876918 \left(\frac{x}{3}\right)^{11} + \epsilon_0(x)
$$

$$
\mathbf{H}_1(x) =
1.909859286 \left(\frac{x}{3}\right)^2 - 1.145914713 \left(\frac{x}{3}\right)^4 \\+ 0.294656958 \left(\frac{x}{3}\right)^6 - 0.042070508 \left(\frac{x}{3}\right)^8 \\+ 0.003785727 \left(\frac{x}{3}\right)^{10} - 0.000207183 \left(\frac{x}{3}\right)^{12} + \epsilon_1(x)
$$

where $$\vert \epsilon_0 \vert < 1.2 \times 10^{-8}$$ and $$\vert \epsilon_1 \vert < 2.5 \times 10^{-9}$$.

For real $$x$$ where $$x \geq 3$$:

$$
\mathbf{H}_0(x)-Y_0(x) = \frac{2}{\pi x}\frac{0.99999906 +  4.77228920\left(\frac{3}{x}\right)^2 + 3.85542044\left(\frac{3}{x}\right)^4 +  0.32303607\left(\frac{3}{x}\right)^6 }{1 +  4.88331068\left(\frac{3}{x}\right)^2 +  4.2895733\left(\frac{3}{x}\right)^4 +  0.52120508\left(\frac{3}{x}\right)^6} + \epsilon_0(x)
$$

$$
\mathbf{H}_1(x)-Y_1(x) = \frac{2}{\pi}\frac{1.00000004 +  3.92205313\left(\frac{3}{x}\right)^2 + 2.64893033\left(\frac{3}{x}\right)^4 +  0.27450895\left(\frac{3}{x}\right)^6 }{1 +  3.81095112\left(\frac{3}{x}\right)^2 +  2.26216956\left(\frac{3}{x}\right)^4 +  0.10885141\left(\frac{3}{x}\right)^6} + \epsilon_1(x)
$$

where $$\vert \epsilon_0 \vert < 8.2 \times 10^{-9}$$ and $$\vert \epsilon_1 \vert < 2.5 \times 10^{-8}$$.

Plotting the absolute error suggests that the error estimations are not overstated, even with the small number of terms used in the approximations:

![Newman abs error]({{ site.url }}/images/2019-11-19-approximating-struve-functions/newman_abs_err.png "Newman abs error")

Although in some places the relative error is somewhat higher than can be achieved using the asymptotic expansion for large values of z:

![Newman rel error]({{ site.url }}/images/2019-11-19-approximating-struve-functions/newman_rel_err.png "Newman rel error")

# Via Quadrature

Struve functions have integral forms. We can directly obtain approximations using quadrature therefore. This is something I'm reluctant to delve into however - the integrands are oscillatory.

# How are other people actually doing it?

As wonderful as consulting outdated papers is, a pragmatic approach would be to follow the others of others. I found several implementations:

* Python: SciPy: A mixed approach. Large z expansion if $$z \geq 0.7v + 12$$, expansion in Bessel functions if $$\vert z \vert < \vert v \vert + 20$$. A power series otherwise.
* Python: mpmath: Using hypergeometric series (power series). I think.
* MATLAB: (Struve function on File Exchange): Chebyshev expansions for small z, and rational approximations for large z - much like a higher accuracy version of Newman.
* C: CEPHES: Hypergeometric series (power series).

# Conclusions (for now)

There appears to be no easy, accurate approximation for the Struve functions. 
To me, the most promising approach seems to be Chebyshev expansions, as used in the MATLAB implementation and Newman. These require effort to produce. A quick implementation probably does best to directly evaluate the power series. This has been my path forward.

# References

* Frank W.J. Olver, Daniel W. Lozier, Ronald F. Boisvert and Charles W. Clark, [NIST Handbook of Mathematical Functions Paperback and CD-ROM](https://isbnsearch.org/isbn/0521140633), Cambridge University Press, 2010

* M. Abramowitz and I. Stegun, [Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables](https://isbnsearch.org/isbn/0-486-61272-4), United States Department of Commerce, National Bureau of Standards, 1964

* Ronald M. Aarts and Augustus J. E. M. Janssen, [Approximation of the Struve function H1 occurring in impedance calculations](https://www.doi.org/10.1121/1.1564019), Acoustical  Society  of America, 2003.

* Ronald M. Aarts and Augustus J. E. M. Janssen, [Efficient approximation of the Struve functions Hn occurring in the calculation of sound radiation quantities](https://www.univie.ac.at/nuhag-php/janssen/data/p227.pdf), 2016.

* J. N. Newman, [Approximations for the  Bessel and Struve Functions, Mathematics of Computation](https://pdfs.semanticscholar.org/2d51/3a7d60e9dfb0b73ea70f2d4d33532973b2f0.pdf), 1984

* R. Zanovello, [Sul calcolo numerico della funzione di Struve Hv⁢  (z)]( 
http://www.seminariomatematico.unito.it/rendiconti/cartaceo/32/Volume%2032%20(1973-74)%20-%2011%20-%20Sul%20calcolo%20numerico%20della%20funzione%20di%20Struve%20Hv(z).pdf), Rend. Sem. Mat. Univ. e Politec. Torino, 1975.


