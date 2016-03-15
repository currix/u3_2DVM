import numpy as np

from matplotlib import pyplot

from survival_probability_basis import read_energy_data

%run histogram_components_LDOS.py

%run energy_basis_elements_so3_chain.py

def sep(xi):
    return (1.-5.*xi)**2/(16.*xi)

# Check states closest to the Separatrix
Nval = 2000
diag_energies_so3 = energy_basis_so3(Nval,np.arange(0,Nval+1,2),0,0.6)

E0 = 366.56652341623544 # N = 2000, L = 0, xi = 0.6

diag_energies_so3 = (diag_energies_so3-E0)/Nval

vvec = (Nval-np.arange(0,Nval+1,2))/2

np.savetxt("fig_ebasis_so3_xi06.dat", np.transpose([vvec, diag_energies_so3]))


wvec = (np.arange(0,Nval+1,2))

np.savetxt("fig_ebasis_so3_xi06.dat", np.transpose([wvec, diag_energies_so3]))

eigenvalues = read_energy_data("../test/eigval_so3_N2000_L0.dat")

diffvec = diag_energies_so3-sep(0.6)

np.abs(diffvec).argmin()
#Out[16]: 578

n=0
cent,bar = hist_components_LDOS("../test/eigvec_so3_N2000_L0.dat", n, "../test/eigval_so3_N2000_L0.dat",N_value=2000, bins=40, OutputCentroid=True, fig=False)
#
diag = hist_components_LDOS("../test/eigvec_so3_N2000_L0.dat", n, "../test/eigval_u2_N2000_L0.dat",N_value=2000, bins=40, OutputCentroid=False, fig=True)
#
np.savetxt("hist_w"+str(n)+"_cent_xi06_N2000_so3.dat",np.transpose(np.array([cent,bar])))
#
np.savetxt("hist_w"+str(n)+"_shape_xi06_N2000_so3.dat",diag)

I

