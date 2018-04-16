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

---



finally, results! 

This is a graph of the results when we keep the mass constant and vary Lambda. On the y-axis we have plot how much the results deviate from the current FRW lensing prediction, in fractional terms. And on the x-axis we have Omega_Lambda, which we are varying. So of course the FRW prediction themselves, the red line, have zero deviations, and the blue points are our numerical results. The green is Rindler and Ishak, and the blue is Kantowski. 

<!-- Towards the right, as Lambda increases, the hole size increases because we're keeping the mass constant while matter density is being reduced.  -->

Our results seem to follow the Kantowski predictions most closely. But there is also a gap between our results and the Kantowski predictions which gets smaller as Lambda increases. If you remember the Kantowski equation for alpha, which you don't, it looks like this:

His analytical calculations were up to some order, and there are higher order terms he neglected, which is this term, the ratio of this term to the lowest order term is around the same order of magnitude as the difference here, and as expected, this term decreases towards the right. So it can probaby explain why our numerical results agree better with the Kantowski predictions at higher Lambda.  

Although the difference here looks quite big, as you can see the y-axis scale is at 10^-7, so this is quite a small effect. 

There are a few factors at play here. 

First, if we keep the mass constant, the size of the hole is related to the matter denity in the universe, and in a flat universe, density decreases as Lambda increases, so the hole gets bigger. So intuitively, alpha increases towards the right because light spends more time in the bigger hole, and all the bending happens inside the hole only. 

Secondly, Lambda also has an effect on how fast the hole boundary expands in static kottler coordinates. So light may have a little more or a little less time in the hole depending on what Lambda is. 

Thirdly, Jacobian at the boundary also has a Lambda dependence. 

The second seems to be a direct effect of how Lambda affects the universe's expansion, but the first one doesn't seem to be a truly direct Lambda effect--it just seems to be an indirect effect of a higher lambda implying a lower matter density, hence causing the hole to be bigger. 

So I did another run for keeping the hole size constant, but in flat space in order to keep the hole size constant we would have to vary the mass. So going towards higher Lambda, for the same hole size, mass would have to decrease. 

So preferably we would want to keep both the hole size and the mass constant and only vary Lambda, but to do that we would have to use a curved universe, which is something that I'm currently working on. 

So the preliminary conclusion is that according to our numerical integrations, the Kantowski results seems to be most correct, GIVEN that the Swiss-Cheese model is an accurate model of lensing in our universe. 

Currently I'm working on doing the FRW calculations in curved space, which would be useful to understanding the true difference that Lambda makes in lensing. Also, numerical errors are currently quite crudely estimated, if possible perhaps it would be good to look into how numerical errors can be properly estimated for the integrator, maybe by varying the step size. 

Lastly, it may be useful to extend the integration to more realistic mass distributions, which would just involve replacing the Kottler hole with another metric and applying the matching conditions again. 

Thank you!

