import numpy as np

from matplotlib import pyplot

from survival_probability_basis import read_energy_data

%run histogram_components_LDOS.py

%run energy_basis_elements_u2_chain.py

def sep(xi):
    return (1.-5.*xi)**2/(16.*xi)

# Check states closest to the Separatrix
Nval = 2000
diag_energies_u2 = energy_basis_u2(Nval,np.arange(0,Nval+1,2),0,0.6)

eigenvalues = read_energy_data("../test/eigval_u2_N2000_L0.dat")
E0 = eigenvalues[0]
diffvec = (diag_energies_u2-E0)/2000-sep(0.6)
np.abs(diffvec).argmin()
#Out[20]: 0
np.abs(diffvec[1:]).argmin()
#Out[21]: 666

#
cent,bar = hist_components_LDOS("../test/eigvec_u2_N2000_L0.dat", 0, "../test/eigval_u2_N2000_L0.dat",N_value=2000, bins=40, OutputCentroid=True, fig=True)
pyplot.show()
#
np.savetxt("hist_n0_cent_xi06_N2000.dat",np.transpose(np.array([cent,bar])))
#
diag = hist_components_LDOS("../test/eigvec_u2_N2000_L0.dat", 0, "../test/eigval_u2_N2000_L0.dat",N_value=2000, bins=40, OutputCentroid=False, fig=True)
#
np.savetxt("hist_n0_shape_xi06_N2000.dat",diag)

I

