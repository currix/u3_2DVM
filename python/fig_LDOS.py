import numpy as np
from matplotlib import pyplot
%run histogram_components_LDOS.py
cent,bar = hist_components_LDOS("../test/test_anharmonicity/eigvec_u2_N600_xi045_a-04.dat", 0, "../test/test_anharmonicity/energ_N600_xi045_a-04.dat",N_value=600, bins=30, fig=True)
pyplot.show()
np.savetxt("hist_n0.dat",np.transpose(np.array([cent,bar])))
# Check states closest to the Separatrix
diffvec = abs(diag_energies_u2 - sep_mean_field_limit(0.45,-0.4,1))
In [138]: diffvec.min()
Out[138]: 0.00010362238774191956

In [139]: diffvec.argmin()
Out[139]: 210  # n = 420

In [140]: diffvec[210] = diffvec[210]+1

In [141]: diffvec.min()
Out[141]: 0.00047790237575279226

In [142]: diffvec.argmin()
Out[142]: 1 # n = 2

In [143]: diffvec[1] = diffvec[1]+1

In [144]: diffvec.min()
Out[144]: 0.00068792956637292457

In [145]: diffvec.argmin()
Out[145]: 0 # n = 0
