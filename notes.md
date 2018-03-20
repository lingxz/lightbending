The 3 distances in lensing: r0, b, R
When we say this is D_L*theta, which are we referring to? Because they are the same to first order but differ in higher orders. 


---


# LTB hole

- Is LTB in comoving or static coordinates?
- How do you tell? Is it determined at the start or do they fall out of the EFEs and then you identify them as static/comoving
- how to determine E and E_r?
- consequently, how to determine R_r?

---

metrics
ds^2 = -dt^2 + a(t)^2 \left [ \frac{dr^2}{1-kr^2} + r^2(d\theta^2 + \sin^2\theta d\phi^2) \right ]
ds^2 = -f(R)dT^2 + \frac{dR^2}{f(R)} + R^2(d\theta^2 + \sin^2 \theta d\phi^2)
f(R) = 1-\frac{2M}{R} - \frac{\Lambda R^2}{3}

geodesic equations
\ddot{x}^{\mu} + \Gamma^{\mu}_{\alpha \beta} \dot{x}^{\alpha} \dot{x}^{\beta} = 0 

frw:
\dot{t} = -\sqrt{\frac{a^2\dot{r}^2}{1-kr^2} + a^2r^2 \dot{\phi}}
a_{,t} = aH_0 \sqrt{\Omega_M/a^3 + \Omega_k/a^2 + \Omega_{\Lambda}}
\ddot{r}  = (1-kr^2)r\dot{\phi}^2 - \frac{k\dot{r}^2}{1-kr^2} - \frac{2a_{,t}}{a}\dot{r}\dot{t}

kottler:
\ddot{R}  = \frac{L_{\text{k}}^2 (R-3M)}{R^4}
\dot{\phi} = \frac{L_k}{R^2}
R_{h,t} = \left ( 1 - \frac{2M}{R_h} - \frac{\Lambda R_h^2}{3} \right ) \sqrt{\frac{2M}{R_h} + \frac{\Lambda R_h^2}{3}}

gravitational lensing equations
D_s \theta_E = D_{LS}\alpha
\alpha_{\text{FRW}} = 4\frac{M}{R} + \frac{15\pi}{4}\left ( \frac{M}{R} \right )^2 + \frac{401}{12}\left ( \frac{M}{R} \right )^3
\alpha_{\text{Ishak}} &= 4\frac{M}{R} + \frac{15\pi}{4}\left ( \frac{M}{R} \right )^2 + \frac{305}{12}\left ( \frac{M}{R} \right )^3 - \frac{\Lambda R r_h}{3}


---

Instead of working directly with Lambda, throughout this paper we work with Omega_Lambda instead.

Selection of the step size is one of the most important concepts in numerical integration of differential equation systems. It is not practical to use constant step size in numerical integration. If the selected step size is large in numerical integration, the computed solution can diverge from the exact solution. And if the chosen step size is small, the calculation time, number of arithmetic operations, and the calculation errors start to increase. So, if the solution is changing rapidly, the step size should be chosen small. Inversely, if the solution is changing slowly, we should choose bigger step size.

Without loss of generality, we consider null geodesics on the plane with $\theta = \pi/2$. 

The Kottler/LTB region has to be mass compensating, meaning the enclosed mass has to equal the mass excised from the FRW background. 

The study of light
propagation in such a background was pioneered by Kantowski

To match
a Szekeres patch to a Friedmann background across a
comoving spherical surface, r = constant, the conditions
are: that the mass inside the junction surface in the Szekeres
patch is equal to the mass that would be inside that
surface in the homogeneous background; that the spatial
curvature at the junction surface is the same in both
the Szekeres and Friedmann models, w

Numerical limitations:

Even higher resolution would be computationally expensive, up to the point where numerical errors in the integration will dominate. 

Nevertheless, a very basic problem arising in cosmology, and which
hinders comparison between theory and observations in other models than the standard Friedmann-Lemaˆıtre-Robertson-Walker (FLRW), stems from the
fact that as soon as we depart from the simple FLRW models largely used by
observational cosmologists, the task of finding null geodesic solutions quickly
becomes an intractable analytical problem. 

The Lemaˆıtre-Tolman-Bondi (LTB) spacetime is the most general spherically
symmetric dust solution of Einstein’s field equations, having the FLRW models
as special sub-cases

# Resources for boundary matching

- https://tel.archives-ouvertes.fr/tel-01235603/document
formula 1.27 on p22

- https://core.ac.uk/download/pdf/25282922.pdf
formula 15 on p4 for extrinsic curvature

# Gravitational lensing for distributed mass

- thin screen approximation

The extent of the mass distribution is very small compared to the distances between source, lens, and observer. Therefore, the mass distribution of the lenns can be treated as if it were an infinitely thin mass sheet perpendicular to the line-of-sight. 

The surface mass is simply obtained by projection. 

The plane of the mass sheet is called the lens plane. 

For an extended mass, it is the PROJECTED mass that counts, not the mass enclosed in that radius. 

Calculating projected mass for NFW: https://arxiv.org/pdf/astro-ph/0002395.pdf
(eq 43, page 7)

When Lambda changes, many things change. What should we keep constant?
