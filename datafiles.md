diff_lambdas_small2: latest, diff z_lens, 100 per Omega_Lambda, but wrong DL
diff_lambdas_small: same z_lens, diff z_source, 100 per Omega_Lambda
diff_lambdas_small3: same as diff_lambdas_small2, but with correct DL
diff_lambdas_small4: same as previous, but the step size in kottler solver is changed to agree with dt. Previously was fixed at 5e-7, oops!!!
diff_lambdas_bigger_redshifts: bigger redshifts, 0.2 to 1 instead of 0.05 to 0.2.
diff_lambdas_bigger_redshifts2: same as above, just more points
diff_lambdas_const_rh: constant rh, varying mass
diff_lambdas_const_rh2: same as above, bigger mass (10^13 solar masses instead of 10^12)
diff_lambdas_const_rh_smaller_redshifts: same as diff_lambdas_const_rh, just smaller redshifts (0.05 - 0.2) instead of 0.2 - 1.


data/lens_z_omlambda_0.0.csv: normal
data/lens_z_omlambda0_2.csv: same as above, but step size in kottler solver is changed to agree with dt
lens_z_omlambda0_2_bigger_redshift.csv: same as above, but redshifts are 0.2 to 1. instead of 0.05 to 0.2
<!-- previously was fixed at 5e-7 -->