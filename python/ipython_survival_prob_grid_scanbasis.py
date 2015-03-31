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
##eigenvectors = "../test/eigvec_so3_xi04_N500_L0.dat"
eigenvectors = "../test/eigvec_u2_xi04_N500_L0.dat"
##eigenvalues = "../test/energy_so3_xi_04_N500_L0.dat"
eigenvalues = "../test/energy_u2_xi_04_N500_L0.dat"
#
## time_max = 16
## time_points = 300
time_max = 0.5
time_points = 200
nvals = np.arange(0, N_value + 2, 2)
#
eigstates=read_eigenstates(eigenvectors)
eigvals= read_energy_data(eigenvalues)
#
time_grid = np.linspace( 0.0, time_max, time_points)
spb = survival_probability_basis_states(eigvals, eigstates, 1, N_value/2 + 1, t_max = time_max, dim_t = time_points)
#
##majorXLocator   = MultipleLocator(5)
##majorXFormatter = FormatStrFormatter('%d')
majorXLocator   = MultipleLocator(0.2)
majorXFormatter = FormatStrFormatter('%3.1f')
#
fig,axes = pyplot.subplots(5,5,sharex=True,sharey=True)
#
pyplot.xlim(0,time_max)
pyplot.ylim(0,1)
#
for index_i in np.arange(0,5):
    for index_j in np.arange(0,5):
        basis_state = 50*index_i + 10*index_j
        #
        axes[index_i,index_j].plot(time_grid,spb[basis_state,:],'b') # U(2)
        ##axes[index_i,index_j].plot(time_grid,spb[basis_state,:],'r') # SO(3)
        ##axes[index_i,index_j].annotate(basis_state, xy=(11, 0.75))
        axes[index_i,index_j].annotate(basis_state, xy=(0.35, 0.75))
        axes[index_i,index_j].xaxis.set_major_locator(majorXLocator)
        axes[index_i,index_j].tick_params(axis='both', which='major', labelsize=9)
