---
layout: post
title: The nasty integral in the Sclavounos unsteady lifting-line theory integro-differential equation
category: posts
draft: false
---

*Edit 03/09/2019: Correction to integral $$I_2$$.*

A foreword: this post is for myself and anyone who is confused at how I obtained the maths used in part of LiftingLineTheory.jl. It also serves as an opportunity for me to check my maths. Consequently, it isn't very exciting. 

---

Sclavounos' unsteady lifting-line theory is attractive because it applies to the frequency range for which finite wing effects are of importance, but doesn't have the complexity (or capability) of that of Guermond and Sellier's 1991 theory.

To obtain a solution, an integro-differential equation must be solved. This gives the user a warm fuzzy feeling since this is similar to the Prandtl method. Sneakier theories (Van Dyke 1964, Guermond 1989, Guermond and Sellier 1991) don't require this step. 

The integro-differential equation is

$$ \Gamma_i(y) - \frac{d_{hn}(y)}{2 \pi i \omega} \int^s_{-s} \frac{d\Gamma_i(\eta)}{dy} K(y - \eta) d\eta = d(y) $$

where the bound circulation $$\Gamma$$ is approximated by a Fourier series $$\Gamma = \sum_{k=1}^N \Gamma_k \sin(k \theta)$$, the spanwise coordinate $$y = s \cos(\theta) \text{ for } \theta \in [0, \pi]$$, $$d(y)$$ is the 2D bound circulation (and can be computed in advance) and $$\omega > 0$$ is angular frequency. $$\Gamma$$ and $$d(y)$$ are complex quantities.

The kernel $$K(y)$$ is given by Sclavounos as 

$$K(y) = \frac{1}{2} \text{sgn}(y) \left[ \frac{e^{-\nu \lvert y \rvert}}{\lvert y \rvert} - i \nu E_1(\nu \lvert y \rvert) + \nu P(\nu \lvert y \rvert ) \right] $$

where $$\nu = \omega / U$$ (where $$U$$ is the free stream velocity, $$E_1(x)$$ is the exponential integral as defined in Olver et al., and 

$$P(y) = \int^\infty_1 dt\;e^{-yt}\left[\frac{\sqrt{t^2 - 1} - t}{t} \right] + i \int^1_0 dt \; e^{-yt}\left[ \frac{\sqrt{1 - t^2} - 1}{t} \right]$$

Sclavounos and his friend Newman are clearly rather sharp to have reached these expressions without the help of a CAS. $$K$$ relates to the spanwise of the problem. For example, for a steady problem it would be $$1/2y$$.

But for the sake of this post, the question is how do we evaluate the integral in the integro-differential equation? Alas, both the $$1/y$$ and the $$E_1$$ terms contain singularities. We want a numerically manageable expression for

$$ I = \int^s_{-s} d\eta \frac{d\Gamma_i(\eta)}{dy} K(y - \eta) $$

Some preliminaries:

$$ \frac{d\Gamma(\eta)}{d\eta} = \frac{d\Gamma(\eta)}{d\theta}\frac{d\theta}{d\eta} $$

so (inserting a minus sign so we can swap the limits)

$$ I = -\int^\pi_0 d\theta \frac{d\Gamma(\eta)}{d\theta} K(y - \eta) $$

where

$$ \frac{d\Gamma(\theta)}{d\theta} = \sum \Gamma_i k \cos(k \theta) $$

This might further be decomposed into three expressions as

$$I = I_1 + I_2 + I_3$$

where

$$ I_1 = -\int^\pi_0 d\theta \frac{d\Gamma(\eta)}{d\theta} \frac{1}{2} \text{sgn}(y-\eta) \frac{e^{-\nu \lvert y-\eta \rvert}}{\lvert y -\eta\rvert} $$

$$ I_2 = \int^\pi_0 d\theta \frac{d\Gamma(\eta)}{d\theta} \frac{1}{2} \text{sgn}(y-\eta)i \nu E_1(\nu \lvert y - \eta \rvert) $$

$$ I_3 = -\int^\pi_0 d\theta \frac{d\Gamma(\eta)}{d\theta} \frac{1}{2} \text{sgn}(y-\eta)\nu P(\nu \lvert y - \eta \rvert ) $$

We'll tackle $$I_1$$ and $$I_3$$ first, since they're a little easier, and then $$I_2$$.

---

Starting with $$I_1$$ we can immediately simplify it

$$ I_1 = -\int^\pi_0 d\theta \frac{d\Gamma(\eta)}{d\theta} \frac{1}{2} \frac{e^{-\nu \lvert y-\eta \rvert}}{y -\eta} $$

and further

$$ I_1 = - \frac{1}{2}\sum \Gamma_i k \left[ \int^\pi_0 d\theta \cos(k \theta) \frac{e^{-\nu \lvert y-\eta \rvert}}{y -\eta} \right]$$

where $$1/(y - \eta)$$ can be identified as the singular part, singular at $$\eta = y = s \cos(\theta_s)$$. So we can evaluate using the singularity subtraction method.

$$ I_1 = - \frac{1}{2}\sum \Gamma_i k \left[ \int^\pi_0 d\theta  \frac{1}{y - \eta} \left( \cos(k\theta) e^{-\nu \lvert y - \eta \rvert} - \cos(k\theta_s)\right) + \cos(k \theta_s) \int^\pi_0 d\theta\;\frac{1}{y - \eta} \right]$$

$$ I_1 = - \frac{1}{2}\sum \Gamma_i k \int^\pi_0 d\theta  \frac{1}{y - \eta} \left( \cos(k\theta) e^{-\nu \lvert y - \eta \rvert} - \cos(k\theta_s)\right)$$

which can now be evaluated numerically. The quadrature should be split at $$\theta_s$$ since due to the discontinuous derivative.

---

Next we'll tackle $$I_3$$ because its easier than $$I_2$$. The good news is that $$P(y)$$ isn't singular! Consequentially we can numerically evaluate 

$$ I_3 = -\frac{\nu}{2} \sum \Gamma_i k\left[ \int^\pi_0 d\theta \cos(k\theta)\text{sgn}(y-\eta) P(\nu \lvert y - \eta \rvert ) \right]$$

again splitting our quadrature about $$\eta = y$$. A minor difficulty is in evaluating $$P(y)$$. Repeating $$P(y)$$,

$$P(y) = \int^\infty_1 dt\;e^{-yt}\left[\frac{\sqrt{t^2 - 1} - t}{t} \right] + i \int^1_0 dt \; e^{-yt}\left[ \frac{\sqrt{1 - t^2} - 1}{t} \right]$$

one finds that the first term seems to require many many quadrature points to obtain a good solution.

We can make the integral more amenable by using integration by parts.

$$
\begin{multline}
\int^\infty_1 dt\;e^{-yt}\left[\frac{\sqrt{t^2 - 1} - t}{t} \right]
 = \\ -e^{-yt}\left(\frac{\pi}{2} -1\right) + \int^\infty_1 dt\;
e^{-yt}\left(\sin^{-1}\left(\frac{1}{\lvert t \rvert }\right)  + \sqrt{t^2 -1} - t\right)
\end{multline}$$

which seems to be easier to evaluate numerically (although, given that I'm working with big floating point numbers, it might be that something funny is happening). I used Gauss-Laguerre quadrature.

---

$$I_2$$ is a challenge.

 
$$ I_2 = \int^\pi_0 d\theta \frac{d\Gamma(\eta)}{d\theta} \frac{1}{2} \text{sgn}(y-\eta)i \nu E_1(\nu \lvert y - \eta \rvert) $$


$$ I_2 = \frac{i \nu}{2} \sum \Gamma_k k \left[\int^\pi_0 d\theta \cos(k \theta)\text{sgn}(y-\eta)E_1(\nu \lvert y - \eta \rvert) \right]$$

Eeek! We can't easily apply singularity subtraction since there isn't an easily identifiable singular bit with a known integral. To get something with a known integral takes a little work. Lets look at something simple and see if we can work backwards. We know that

$$ \frac{dE_1(z)}{dz} = - \frac{e^{-z}}{z}$$

Working from this, a simple integral is

$$  \int dx\;E_1(x) = x E_1(x) - e^{-x} + C $$

$$  \int dx\;E_1(-x) = x E_1(-x) + e^{x} + C $$

where I've used integration by parts.
Extending this concept a little in the Cauchy sense with a singularity at $$x=0$$, where $$a<0<b$$ :

$$\int^b_a dx \; \text{sgn}(x)E_1(\lvert x \rvert) = 
\lim_{\varepsilon\rightarrow0^+} \left( \left[-xE_1(-x) - e^x \right]^{-\varepsilon}_a + \left[xE_1(x) -  e^{-x} \right]^b_{\varepsilon} \right) $$

$$\int^b_a dx \; \text{sgn}(x)E_1(\lvert x \rvert) = 
\lim_{\varepsilon\rightarrow0^+} \left(
\left[ \varepsilon E_1(\varepsilon) - e^{-\varepsilon} +
 a E_1(-a) + e^{a} \right] + 
\left[ bE_1(b) - e^{-b}  
-\varepsilon E_1(\varepsilon) + e^{-\varepsilon} \right]
 \right) $$

Since $$\lim_{\varepsilon\rightarrow0^+}(\varepsilon E_1(\varepsilon)) = 0$$, this evaluates as 

$$
\int^b_a dx \; \text{sgn}(x) E_1(\lvert x \rvert) = aE_1(-a) + e^a + bE_1(b) - e^{-b}
$$

We can substitute some of the complexity of the original kernel back in as $$x = \nu (y - \eta)$$, swapping the integration variable to $$\theta$$ again.
Consequently,

$$
x = \nu ( y - s \cos(\theta) )
$$

$$
dx = \nu s \sin(\theta)\; d\theta
$$

$$
\begin{multline}
\int^b_a dx \; \text{sgn}(x)E_1(\lvert x \rvert) = \\
\nu s \int^\pi_0 d\theta \; \text{sgn}( y - s \cos(\theta))E_1(\nu \lvert y - s \cos(\theta) \rvert) \sin(\theta)
\end{multline}
$$

where, to satisfy the limits,
$$b = \nu ( y - s \cos(\pi)) = \nu(y + s)$$ and $$a = \nu(y - s)$$.
So,

$$
\begin{multline}
\nu s \int^\pi_0 d\theta \; \text{sgn}( y - s \cos(\theta))E_1(\nu \lvert y - s \cos(\theta) \rvert) \sin(\theta) = \\
\nu(y - s)E_1(-\nu(y - s)) + e^{\nu(y - s)} +
\nu(y + s)E_1(\nu(y + s)) - e^{-\nu(y + s)}
\end{multline}
$$

Now we're getting somewhere. For brevity,
if we define

$$F_s(\theta) = \nu s \text{ sgn}( y - s \cos(\theta))E_1(\nu \lvert y - s \cos(\theta) \rvert) \sin(\theta) $$

and 

$$G(\theta) = \frac{\cos(k \theta)}{\nu s \sin(\theta)}$$

then 

$$ \int^\pi_0 d\theta \; G(\theta) F_s(\theta) = \cos(k \theta) \text{ sgn}(y - \eta) E_1(\nu \lvert y - \eta \rvert ) $$

Using the singularity subtraction method,

$$ \int^\pi_0 d\theta \; G(\theta) F_s(\theta) = \int^\pi_0 d\theta \; F_s(\theta) \Big(G(\theta) - G(\theta_s) \Big) + G(\theta_s) \int^\pi_0 d\theta \; F_s(\theta)  $$

where $$\theta_s = \cos^{-1}(y / s)$$. The regularised integral can be evaluated numerically. 

So finally we obtain an admittedly rather cumbersome expression,

$$
\begin{multline}
I_2 = \frac{i \nu}{2} \sum \Gamma_k k \Bigg[\\
\int^\pi_{\theta_s} d\theta \; 
E_1(\nu ( y - s \cos(\theta) )) \sin(\theta) \left(
\frac{\cos(k \theta)}{\sin(\theta)} - \frac{\cos(k \theta_s)}{\sin(\theta_s)}\right) - \\
\int^{\theta_s}_0 d\theta \; 
E_1(\nu ( s \cos(\theta) - y )) \sin(\theta) \left(
\frac{\cos(k \theta)}{\sin(\theta)} - \frac{\cos(k \theta_s)}{\sin(\theta_s)}\right) + \\
\frac{\cos(k \theta_s)}{\nu s \sin(\theta_s)}
\left(
\nu(y - s)E_1(-\nu(y - s)) + e^{\nu(y - s)} +
\nu(y + s)E_1(\nu(y + s)) - e^{-\nu(y + s)}\right)\\
\Bigg]
\end{multline}
$$

---

I do wonder if there is a better way of doing  $$I_2$$ since the current method is a little awkward - sufficiently awkward that I shouldn't be surprised if I've made a mistake somewhere. However, this is for better or worse the method currently in use by my implementation in LiftingLineTheory.jl.

---

Recommended reading:

* P. D. Sclavounos, [An unsteady lifting-line theory](https://www.doi.org/10.1007/bf00127464), Journal of Engineering Mathematics, 1987
* Milton Van Dyke,[Lifting-line theory as a singular-perturbation problem](https://www.doi.org/10.1016/0021-8928(64)90134-0), Journal of Applied Mathematics and Mechanics, 1964 
* Jean-Luc Guermond and Antoine Sellier, [A unified unsteady lifting-line theory](https://www.doi.org/10.1017/s0022112091003099), Journal of Fluid Mechanics, 1991 
* Frank W.J. Olver, Daniel W. Lozier, Ronald F. Boisvert and Charles W. Clark, [NIST Handbook of Mathematical Functions Paperback and CD-ROM](https://isbnsearch.org/isbn/0521140633), Cambridge University Press, 2010