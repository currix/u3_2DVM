# IPython log file
import numpy as np
from matplotlib import pyplot
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#
from survival_probability_basis import survival_probability_basis_states, read_energy_data, read_eigenstates, energy_basis_I, energy_basis_II
#
N_value = 500
l_value = 0
xi_value = 0.4
eigenvectors = "../test/eigvec_so3_xi04_N500_L0.dat"
##eigenvectors = "../test/eigvec_u2_xi04_N500_L0.dat"
eigenvalues = "../test/energy_so3_xi_04_N500_L0.dat"
##eigenvalues = "../test/energy_u2_xi_04_N500_L0.dat"
#
time_max = 15
time_points = 300
#
eigstates=read_eigenstates(eigenvectors)
eigvals= read_energy_data(eigenvalues)
nvals = np.arange(0, N_value + 2, 2)
#
energy_bas_I = energy_basis_I(N_value, nvals, l_value, xi_value)
energy_bas_II = energy_basis_II(N_value, nvals[::-1], l_value, xi_value)
#
time_grid = np.linspace( 0.0, time_max, time_points)
spb = survival_probability_basis_states(eigvals, eigstates, 1, N_value/2 + 1, t_max = time_max, dim_t = time_points)
#
majorLocator   = MultipleLocator(5)
majorFormatter = FormatStrFormatter('%d')
minorLocator   = MultipleLocator(1)
#
## for plot in (nvals/2)[1::10]:
for plot in range(52,61):
    #
##    panel = pyplot.subplot(5, 5, plot/10+1)
    panel = pyplot.subplot(3, 3, plot-51)
    #
    pyplot.ylim(0,0.4)
    pyplot.annotate(plot, xy=(8.5, 0.3))
    #
    pyplot.plot(time_grid,spb[plot,:],'r')
    #
    panel.xaxis.set_major_locator(majorLocator)
    panel.xaxis.set_major_formatter(majorFormatter)
    panel.xaxis.set_minor_locator(minorLocator)
    panel.tick_params(axis='both', which='major', labelsize=9)
