diff_lambdas_small2: latest, diff z_lens, 100 per Omega_Lambda, but wrong DL
diff_lambdas_small: same z_lens, diff z_source, 100 per Omega_Lambda
diff_lambdas_small3: same as diff_lambdas_small2, but with correct DL
diff_lambdas_small4: same as previous, but the step size in kottler solver is changed to agree with dt. Previously was fixed at 5e-7, oops!!!
diff_lambdas_bigger_redshifts: bigger redshifts, 0.2 to 1 instead of 0.05 to 0.2.
diff_lambdas_bigger_redshifts2: same as above, just more points

diff_lambdas_const_rh: constant rh, varying mass
diff_lambdas_const_rh2: same as above, bigger mass (10^13 solar masses instead of 10^12)
diff_lambdas_const_rh_smaller_redshifts: same as diff_lambdas_const_rh, just smaller redshifts (0.05 - 0.2) instead of 0.2 - 1.
diff_lambdas_const_rh_smaller_step: same as diff_lambdas_const_rh (10^12 solar masses, redshifts 0.2-1), but with smaller step size 1e-7
dopri5_diff_lambdas_const_rh_smaller_step [INTERRUPTED]: same as above, but using dopri5 instead of vode, redshifts change to 0.1 to 0.5 instead of 0.2-1. 
diff_lambdas_const_rh3 [INTERRUPTED]: same as diff_lambdas_const_rh, but bigger mass, higher theta (10^13 solar masses, 5e-5 rad), redshifts 0.1-0.5, step size 1e-7

diff_lambdas_ltb_cheese: static general mass distribution inside the hole

data/lens_z_omlambda_0.0.csv: normal
data/lens_z_omlambda0_2.csv: same as above, but step size in kottler solver is changed to agree with dt
lens_z_omlambda0_2_bigger_redshift.csv: same as above, but redshifts are 0.2 to 1. instead of 0.05 to 0.2
<!-- previously was fixed at 5e-7 -->


diff_lambdas_redshifts: from sw_lambda3.py. different redshifts, different lambda, same mass, fixed start theta. 
diff_lambdas_masses: from sw_lambda4.py. different masses, same redshift, different lambda, fixed start theta.