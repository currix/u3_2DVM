def read_energy_data(filename):
    #
    import numpy as np
    #
    # Read fortran output
    fortran_output = np.loadtxt(filename)
    #
    # Select secod column, eigenvalues
    eigenvalues = fortran_output[:,1]
    #
    return eigenvalues
#########################################################################
#########################################################################
def survival_probability(Eigenvalues, alpha_values, index = 0, t_max = 3, dim_t = 100):
    '''Compute the survival probability of the eigenstate superposition phi_n. See notes... (TODO complete docstring)'''
    #
    import numpy as np
    #
    time_grid = np.linspace( 0.0, t_max, dim_t)
    #
    if (index == 0):
        num_curves = (len(Eigenvalues)-1)/3
        #
        result = np.zeros([num_curves,dim_t])
        #
        for curve in range(1,num_curves+1):
            #
            index_n = 3*(curve - 1)
            #
            result[curve-1,:] = compute_surv_prob(Eigenvalues[index_n + 1:index_n + 4], alpha_values, time_grid) 
            #
            #
    return result
#
def compute_surv_prob(Energies, alpha_values, time_grid):
    #
    import numpy as np
    #
    surv_prob = np.zeros(len(time_grid))
    for index in range(3):
        surv_prob = surv_prob + alpha_values[index]**2*np.exp(Energies[index]*1j*time_grid)
    #
    return abs(surv_prob)**2 
#########################################################################
#########################################################################
def plot_surv_prob(surv_prob_output, t_max = 3, dim_t = 100, esqpt_index = -1):
    #
    import numpy as np
    #
    from matplotlib import pyplot
    #
    time_grid = np.linspace( 0.0, t_max, dim_t)
    #
    curve_number = surv_prob_output.shape[0]
    #
    pyplot.ylabel('Survival Probability')
    pyplot.xlabel('time')
    #
    for index in range(curve_number):
        if (index < esqpt_index):
            line_style = "r--"
        elif (index == esqpt_index):
            line_style = "b-"
        else:
            line_style = "g-."
        pyplot.plot(time_grid, surv_prob_output[index,:], line_style)
#########################################################################
#########################################################################

                    
