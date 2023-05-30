## simulate_and_gds.py

#Basic idea:
# 1. Take a cavity parameter array cp, simulate the cavity to get freq, vmode, qe, qi
# 2. Parse into phidl using chang's code (SiC_gds_gen_v5.py, SiC_Overcoupled_gds.ipynb, SiC_cavity_array_code.ipynb)
# 3. Spit out gds of full cavity and a line with parameter values, frequency, qe, qi
