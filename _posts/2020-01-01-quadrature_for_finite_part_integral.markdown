---
layout: post
title: Quadrature for finite part integrals
category: posts
draft: false
---

Recently I've had to work with Hadamard finite part integrals as a result of looking at Guermond and Sellier's 1991 unsteady lifting-line theory. A finite part integral typically looks a little like this:

$$I = \text{F.P.}\int_a^b \frac{f(x)}{(x-s)^2} \;dx$$

where $$a < s < b$$ and $$\text{F.P.}$$ indicates that its a finite part integral (although in the typically useful fashion of mathematical notation, there are a few ways of writing it).
What is immediately apparent is that this integral is divergent so long as $$f(x)$$ is continuous. Fortunately - as the title implies - we only take the finite part.

To me, taking only the finite part is a little confusing - isn't there an infinite selection of finite parts that we could select? I think then, that the following 
form is a little clearer:

$$I = \lim_{\varepsilon\rightarrow0}\left[ \left(\int^b_{s+\varepsilon} + \int^{s-\varepsilon}_a \right) \frac{f(x)}{(x-s)^2} \;dx - \frac{f(s+\varepsilon) + f(s-\varepsilon)}{\varepsilon} \right]$$

Initially I set about integrating and taking limits then. Which is all very time consuming, and not terribly adaptable. So I was pleased to find a numerical method, which I'm posting about here because its so darned neat.

# Paget's finite part quadrature

An alternate way of looking at a finite part integral is to think of it as the derivative of a Cauchy principal value integral:

$$I = \frac{d}{ds} \text{C.P.V.}\int^b_a \frac{f(x)}{(w-s)^2} \; dx$$

Integrating C.P.V. integrals numerically is well understood, with singularity subtraction probably being the most common method.
Lets make the C.P.V. form of $$I$$ a little more convoluted by adding a weight function, $$w(x)$$.

$$I = \frac{d}{ds} \text{C.P.V.}\int^b_a \frac{w(x)g(x)}{(w-s)^2} \; dx$$

The weight function might be familiar from Gaussian quadrature. For an integral $$\int_a^b w(x)h(x) dx$$ the polynomial orthogonal about $$w(x)$$ can be used to compute a quadrature that is exact so long as $$h(x)$$ is a polynomial. The $$n$$ points and weights of this polynomial will be referred to as $$x_i$$ and $$\mu_i$$ respectively. Hence $$\int_a^b w(x)h(x) dx = \sum_{i=1}^n \mu_i h(x_i) $$. If we apply this quadrature to $$I$$ along with singularity subtraction for a singularity at $$s$$, we obtain

$$I = \frac{d}{ds}\left[\sum_{i=1}^n \mu_i\frac{g(x_i) - g(s)}{x_i-s} + q_0(s)g(s) \right] $$

where $$q_0(s) = \int^b_a \frac{w(x)}{x-s} dx$$.

Presented so, its clear how a quadrature for a finite part integral might be obtained. Differentiating:

$$I = \sum_{i=1}^n \mu_i\left[\frac{g(x_i) - g(s)}{(x_i-s)^2} - \frac{g'(s)}{(x_i-s)}\right] + q'_0(s)g(s) + q_0(s)g'(s) $$

where $$\bullet'$$ indicates the derivative with respect to $$s$$. 

Its worth noting that the usual rules on $$g(x)$$ apply - for the above to be exact, $$g(x)$$ should be a polynomial of degree $$2n-1$$, but good results can be obtained even if it isn't. Being polynomial like is important however, and obtaining good accuracy is challenging if the $$g(x)$$ does not have a bounded derivative.

For example, in my case I had $$f(x) = \sin(m \cos^{-1}(x))$$, which has an unbounded first derivative at $$x = \pm 1$$. Consequently Gauss-Legendre quadrature ($$w(x) =1$$) did not produce good results. 
Better results were obtained using a type 2 Gauss-Chebyshev quadrature, meaning $$w(x)=\sqrt{1-x^2}$$ and $$g(x) = \sin(m \cos^{-1}(x)) / \sqrt{1-x^2}$$.

# References

* D.F. Paget, [The numerical evaluation of Hadamard finite-part integrals](https://doi.org/10.1007/BF01395957), Numer. Math., 1981

