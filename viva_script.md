1696 words

So my project is about investigating the role that the cosmological constant Lambda plays in gravitational lensing. 

The key question we want to answer is: does the cosmological constant directly affect gravitational lensing? 

Well there are several reasons why we want to know this. 

The cosmological constant is an important part of our current model of the universe, and gravitational lensing is one of the key observational tests of General Relativity, so we definitely want to know how these two interact with each other. 

And also if we do find an effect, which is generally agreed to be second order or higher, if any, it can become important for future cosmology measurements, although we might not be able to detect it at the moment since Lambda is small.

So an overview of the literature. There are two key papers which basically separates the research into two camps. The first is the conventional view by Islam, which says, no, Lambda does not directly affect gravitational lensing. 

Then about a decade ago Rindler and Ishak published a paper that says, yes, it does. 

And this started the debate and several papers were published after that supporting either camp. 

So this is a picture of the problem--gravitational lensing, and the lens eqution is comes from simple trigonometry, which says that the Einstein angle that we observe depend on the distances to the lens and the source, and the bending angle alpha. 

So the current theory has an expression for alpha which only depends on the mass. These are angular diameter distances, and in the Friedmann Robertson Walker metric they already depend on Lambda, but that is already taken into account. 

So when I say we want to know whether Lambda affects lensing directly, I mean whether there is an extra term we need to take care of. And Rindler and Ishak thinks there is an extra term. So these are the things we want to check our results against. 

Our approach is to do a numerical integration of the paths of the light rays in a Swiss-Cheese model, which I will explain what it is in moment. 

Most of the work done previously has been analytic, but there has been a few numerical works, such as Schucker, who took a partially numerical approach, and Aghili who did numerical integrations in a different metric. 

So in the Swiss Cheese model, we have the cheese outside, which is homogeneous and expanding, described by the Friedmann Robertson Walker metric, and we take a comoving sphere, and collapse all that mass inside to a point, and we have a vacuum that includes a cosmological constant, and the point mass at the centre, and this is described by a Kottler metric. So this hole will be expanding in static coordinates. So the Swiss-Cheese is just two metrics stitched together the boundary with some matching conditions. 

So we propagate the light backwards through the swiss cheese, from the observer back to the source. 

The light travels through the FRW region, then through the Kottler region, then back out into the FRW region. 

To find the light path in each region, we need to first write down the metric, then calculate the christoffels, and get the differential equations using the geodesic equations and, because we're dealing with light, the null condition. Then once we have thoe differential equations, we need to solve that, either analytically or numerically. 

First, in the cheese, the homogeneous part, geometry is described by a FRW metric, from which we get a bunch of differential equations. 

And in the hole, we have the Kottler metric, which is just the Schwarzschild metric extended to include a cosmological constant. 

And in the middle, we have to glue the two metrics together on a hyper surface. The Israel junction conditions give us two criteria: that the induced metric and the extrinsic curvature must match on the boundary. 

Applying these junction conditions, we get a few criteria. The first is the quite intuitive criteria that the value of the point mass inside the hole must be the same as the original mass inside the comoving sphere. 

They also tell us how the boundary of the hole, R_h, expands in the static Kottler coordinates, and also the Jacobian for transforming the velocities from FRW into Kottler coordinates at the boundary. 

So with this we have the full picture. We start off the light ray with a fixed angle from the observer, then propagate the light rays until it reaches the boundary of the hole. Then at the boundary, we convert the velocities from FRW coordinates into Kottler coordinates, and continue propagating the light rays inside the hole. 

Inside the hole, the light ray is moving, and the boundary of the hole is also moving with a different differential equation. At every point we check if it has exited the hole. When it has, we stop the integration and convert back into the FRW, and then we record the coordinate at which it crosses the axis. 

finally, results! This is a graph of the results when we keep the size of the hole constant and vary Lambda. So on the y-axis is how much the results deviate from the current FRW lensing prediction, in fractional terms, and on the x-axis it is plotted against Omega_Lambda. So of course the FRW predictions themselves, the red line, have zero deviations, and the blue point are our numerical results. The green one is the prediction by the Rindler and Ishak, the yes camp people. In this graph, the size of the hole is kept constant, and so the mass decreases as Lambda increases, since the matter density decreases as Lambda goes higher. The error bars are obtained from repeating the integration for different lens redshift.  

This is another run, where we kept the central mass constant. This means that as we increase Lambda the size of the hole increases. 

Our results do not favour the Rindler and Ishak predictions. There is also a downward trend that I am still looking into, I think from looking at the Kantowski paper that did an analytical estimation of the effect, it might be due to the increase in size of the holes, so on this side (towards the right) the light ray spends more time in the hole and hence it's bent more. 

So that is something I am currently working on. 

So that's the first bullet point, comparing with the analytical estimation of the effect in the Swiss-Cheese universe. 

Also, currently the calculations done in the FRW only apply for flat spce, so it would be good to extend that to curved space. Numerically it's not difficult, but since the light ray spends most of its time in the FRW universe, that is going to introduce a lot more numerical errors. 

Also lastly, gravitational lensing could be investigated in a more realistic scenario, for example a NFW mass distribution instead of a point mass. This would just involve replacing the Kottler hole with another LTB metric for a general mass distribution, and applying the matching conditions again. 

Thank you!

