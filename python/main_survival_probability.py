import sys
import numpy as np
from matplotlib import pyplot
from survival_probability import survival_probability, plot_surv_prob, read_energy_data
#
#  python main_survival_probability energy_data_file
#
def main():
    #
    # Extract energy filename from the arguments
    energy_file = sys.argv[1]
    eigval = read_energy_data(energy_file)
    #
    # Define alpha
    alpha = np.array([1.0/2.0,1/np.sqrt(2),1./2.0])
    #
    # Plot survival probability
    sp_results = survival_probability(eigval,alpha)
    #
    # Define Figure
    pyplot.figure(figsize=(20.0,6.0))
    # Plot Figure
    plot_surv_prob(sp_results, esqpt_index = 22)
    # Save Figure
    pyplot.savefig('foo.png', bbox_inches='tight')
    # Display Figure
    pyplot.show()
###################################################################
main()
