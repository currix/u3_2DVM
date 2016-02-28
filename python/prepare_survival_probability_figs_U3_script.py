# Read libraries and functions
import numpy as np
from energy_basis_elements_u2_chain import energy_basis_u2
from energy_basis_elements_so3_chain import energy_basis_so3
from survival_probability_basis import read_energy_data, read_eigenstates, survival_probability_basis_states

#########################################
# Definitions ## Edit only this section #
#########################################

# Define Sytem Size, Angular Momentum, Control Parameter, and Basis
Nval = 1000
Lval = 0
xi = 0.6
basis = "u2" # Options: "u2" or "so3"

## Linear time grid data
t_max = 10;dim_t = 10000
time_grid = np.linspace( 0.0, t_max, dim_t)
## Variable length grids
# t_max = 10000.;dim_t = 10000
# time_grid = np.zeros(dim_t+1)
# for index in np.arange(0,dim_t+1):
#     time_grid[index] = (t_max/dim_t**3)*index**3

# Basis states indexes (start in 1)  -1 --> Closest to separatrix 
basis_states_indexes = [1,2,-1]
###########################################
#  End of Defs  ## Edit only this section #
###########################################
#
print "N = ", Nval, ", L = ", Lval, ", xi = ", xi
## Read eigenvalues (absolute eigenvalues) and eigenvectors
Eigenvalues = read_energy_data("../test/eigval_"+basis+"_N"+str(Nval)+"_L"+str(Lval)+".dat")
Eigenstates = read_eigenstates("../test/eigvec_"+basis+"_N"+str(Nval)+"_L"+str(Lval)+".dat")
E0 = Eigenvalues[0]
Eigenvalues = Eigenvalues - E0 # Excitation energies
# Computing survival probability
print "Computing survival probability ... "
#
sp=np.zeros((time_grid.shape[0],len(basis_states_indexes)+1),order="F") # +1 for the abscyssa values
#
for index in basis_states_indexes:
    if (index == -1):
        # Find basis state close to separatrix
        if basis is "u2":
            basisen = energy_basis_u2(Nval,np.arange(Lval,Nval+1,2), Lval, xi)
        else:
            basisen = energy_basis_so3(Nval,np.arange(Lval,Nval+1,2), Lval, xi)
        #
        ecrit = (1-5*xi)**2/(16*xi)
        # Note that basisen is not normalized
        diffen = (basisen-E0)/Nval-ecrit
        ##
        crit_index = np.abs(diffen).argmin() + 1
        if (crit_index <=1):
            crit_index = np.abs(diffen[2:]).argmin()+3
        print "crit_index = ", crit_index 
        sp[:,index] = survival_probability_basis_states(Eigenvalues,Eigenstates,crit_index,crit_index,time_grid=time_grid)[0]
    else:
        sp[:,index] = survival_probability_basis_states(Eigenvalues,Eigenstates,index,index,time_grid=time_grid)[0]
##
# Add time grid
sp[:,0] = time_grid
print "\t [Done]"

# Save data
print "Saving survival probability ... "
np.savetxt("script_fig_U3_SP_N"+str(Nval)+"_L"+str(Lval)+".dat", sp)
print "\t [Done]"

