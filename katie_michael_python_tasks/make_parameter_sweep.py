# make_parameter_sweep.py
# Using simulate_and_gds, make an array of devices and output a full gds that is ready for writing + a parameter .csv
# This .csv will contain common values at the top, and for each device will have a row with:
# [index, swept cavity parameter values, freq_sim, vmode_sim, qe_sim, qi_sim, maybe some sensitivity stuff]
#
# 1. 'Brute force': vary a large range over many parameters and just make a lot of devices, doing good design <->
# feedback for future rounds
#
# 2. Sensitivity vs. fab uncertainty: For a given 'optimal' design, calculate the
# sensitivity to a range of parameters. If we have sigma_x = delta, s.t. Q_(hx + delta) = Q_thresh,
# then we essentially want to make a parameter sweep s.t. we cover (hx-delta, hx + delta) even with our fabrication
# uncertainty (if I design hx, what is the range of hx that I wouldn't be surprised that I get back)
#
# This idea just came to me in the shower, so we need to flesh it out a bit, but I think that it's the way to go
#
# Requires some knowledge of fabrication uncertainty, maybe what we can do is establish these numbers in fabrication.
# Make a gds with 10-15 of the same device (param_nom). We measure via SEM (param_meas) and can calculate the offset
# from nominal parameters param_off = mean(param_nom - param_meas) and the variance of device parameters param_var =
# var(param_off - param_measure). We can use param_off as a simple fabrication offset, but we should incorporate
# param_mean into approach 2.
