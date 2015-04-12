def plot_components(eigvec_file, eigval_file, bas_state_index, majorX = 50, majorY = 0.02, col="b-o"):
    ###
    '''Plot squared components in all eigenvectors of a bas_state_index basis element (min bas_state_index = 0) as a function of the state energy.''' 
    ###
    import numpy as np
    from matplotlib import pyplot
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    #
    from survival_probability_basis import read_energy_data, read_eigenstates
    #
    #
    eigstates=read_eigenstates(eigvec_file)
    eigvals= read_energy_data(eigval_file)
    #
    #
    majorXLocator   = MultipleLocator(majorX)
    majorYLocator   = MultipleLocator(majorY)
    #
    ##    for bas_state in range(0,4):
    components = eigstates[bas_state_index,:]  ## Components of a given basis element in all eigenstates
    ##
    panel = pyplot.subplot(1, 1, 1)
    panel.tick_params(axis='both', which='major', labelsize=14)
    pyplot.plot(eigvals,components**2,col)
    ##
    panel.xaxis.set_major_locator(majorXLocator)
    panel.yaxis.set_major_locator(majorYLocator)
    ##
    ##
    ##    pyplot.annotate(bas_state, xy=(30, 0.07))
    pyplot.xlabel("Energy")
    pyplot.ylabel("Squared Eigenvector Components")
    panel.xaxis.label.set_fontsize(18)
    panel.yaxis.label.set_fontsize(18)
