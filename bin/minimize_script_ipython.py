import numpy as np
from lmfit import minimize, Parameters, fit_report
from subprocess import Popen, PIPE
#
%run residuals_U3.py
#
cd ../test/test_chi2
# Parameters
#
input_filename = "HNC_input_file"
bin_path = "~/PR_Fortran/U3/PureBending/triat_U3_3.0/bin"
BENT = ".F."
exp_data_file="exp_HNC_Danielle.dat"
N_val=60
LMAX=7
VMAX=7
EMINL=".F."
#
# Hamiltonian Parameters
U3H_FIT_param=Parameters()
U3H_FIT_param.add('P11',value=2378.0, vary=True)
U3H_FIT_param.add('P21',value=-38.0,vary=True)
U3H_FIT_param.add('P22',value=19.0, vary=True)
U3H_FIT_param.add('P23',value=-9.5, vary=True)
U3H_FIT_param.add('P31',value=0.0,vary=False)
U3H_FIT_param.add('P32',value=0.01,vary=True)
U3H_FIT_param.add('P33',value=0.0,vary=False)
U3H_FIT_param.add('P41',value=0.0,vary=False)
U3H_FIT_param.add('P42',value=0.0,vary=False)
U3H_FIT_param.add('P43',value=0.0,vary=False)
U3H_FIT_param.add('P44',value=0.0,vary=False)
U3H_FIT_param.add('P45',value=0.0,vary=False)
U3H_FIT_param.add('P46',value=0.0,vary=False)
U3H_FIT_param.add('P47',value=0.0,vary=False)

out_0 = minimize(residuals_U3, U3H_FIT_param, args=(input_filename, BENT, exp_data_file, N_val, LMAX, VMAX, EMINL))

print fit_report(out_0)
