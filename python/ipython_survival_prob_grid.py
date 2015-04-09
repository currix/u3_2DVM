# IPython log file
import numpy as np
from matplotlib import pyplot
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#
from survival_probability_basis import survival_probability_basis_states, read_energy_data, read_eigenstates
#
N_value = 500
l_value = 0
xi_value = 0.4
##eigenvectors = "../test/eigvec_so3_xi04_N500_L0.dat"
eigenvectors = "../test/eigvec_u2_xi04_N500_L0.dat"
##
eigenvalues = "../test/energy_u2_xi_04_N500_L0.dat"
#
time_max = 1.0
time_points = 200
##time_max = 16
##time_points = 300
#
eigstates=read_eigenstates(eigenvectors)
eigvals= read_energy_data(eigenvalues)
#
time_grid = np.linspace( 0.0, time_max, time_points)
spb = survival_probability_basis_states(eigvals, eigstates, 1, N_value/2 + 1, t_max = time_max, dim_t = time_points)
#
majorXLocator   = MultipleLocator(0.2)
majorXFormatter = FormatStrFormatter('%3.1f')
##majorYLocator   = MultipleLocator(0.2)
##majorYFormatter = FormatStrFormatter('%3.1f')
##
##minorLocator   = MultipleLocator(4)
#
for basis_state in range(0,6): # 
    #
    panel = pyplot.subplot(3, 2, basis_state+1)  ## 
    #
    pyplot.ylim(0,1)
    ##
    pyplot.annotate(basis_state, xy=(0.1, 0.2))
    #
##    pyplot.plot(time_grid,spb[basis_state,:],'r') # SO(3)
    pyplot.plot(time_grid,spb[basis_state,:],'b') # U(2)
    #
    panel.xaxis.set_major_locator(majorXLocator)
    panel.xaxis.set_major_formatter(majorXFormatter)
##    panel.yaxis.set_major_locator(majorYLocator)
##    panel.yaxis.set_major_formatter(majorYFormatter)
##    panel.xaxis.set_minor_locator(minorLocator)
    panel.tick_params(axis='both', which='major', labelsize=9)
