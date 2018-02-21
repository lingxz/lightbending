# Abstract

# Introduction

Describe the problem

Overview of the literature, starting with Rindler and Ishak's paper. Point out other work done on the Swiss Cheese (Kantowski, analytical) and also other numerical work that used a different metric (Aghili, Effect of Accelerated Global Expansion on Bending of Light, used a McVittie metric). 

Talk about the metric, Einstein field equations, and how light path follows the null geodesic. 

Introduce the structure of the rest of the report. 

# Gravitational lensing formalism

## Lensing in the Schwarzschild metric

Derivation of the Schwarzschild lensing angle alpha. 

## Observables

[Diagram of lensing]

Derivation of Einstein angle. 

Make a note that angular diameter distances already depend on Lambda, the question is whether there is an additional dependence on Lambda that we need to take care of (i.e. do we need to modify the existing formula for lensing). 

# Description of the Swiss Cheese model

Describe the Swiss Cheese model. Note that the mass has limited influence, i.e. the bending stops after exiting the hole, back into the cheese. 

Reasons for using this model: is an exact solution of Einstein equations, influence of the central mass confined to the size of the hole

## Spacetime patches

FRW geometry: Describe FRW metric, for general curved space
Kottler geometry: Kottler metric

## Matching conditions

Matching of induced metric and extrinsic curvature
Conversion of velocity between FRW and Kottler spacetime

[Diagram]

## Light propagation

Friedmann equations, null condition, geodesic equations

## Calculation of Einstein angle

How the Einstein angle is calculated from the numerical angular diameter distances

## Hole with a generalised mass distribution

Introduce LTB metric
similar junction conditions at the boundary due to birchoff's theorem

### Density profile

NFW mass profile and lensing mass calculated from thin lens approximation

# Results and discussion

## Numerical integration for a Kottler Swiss Cheese model

Back propagation. 

To check that it's working, it is compared with the Schwarzschild case. 

selection of step size by varying step size for known Schwarzschild case and see how well it agrees with expected result

Choice of numerical values of M (10^12 solar masses) and z_lens

talk about lensed vs unlensed distances: calculate magnitude of fractional uncertainty from Luke Butcher's paper

## Results for a LTB swiss cheese

plots for P(r), rho(r), mass(r) for the NFW density distribution
Results are checked against the Schwarzschild / Kottler case when the mass becomes a point mass at the centre. 


# Conclusion

Summarize what's above
Further work