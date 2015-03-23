def plot_components(eigvec_file, eigstate_index, majorX = 50, majorY = 0.02, col="b-o"):
    ###
    '''Plot the squared components of the eigstate_index-th eigenvector (min value = 0).'''
    #
    import numpy as np
    from matplotlib import pyplot
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    #
    from survival_probability_basis import read_eigenstates
    #
    #
    eigstates=read_eigenstates(eigvec_file)
    #
    majorXLocator   = MultipleLocator(majorX)
    majorYLocator   = MultipleLocator(majorY)
    majorXFormatter = FormatStrFormatter('%d')
    #
    components = eigstates[:,eigstate_index]    ## Eigenstate's components
    ##
    panel = pyplot.subplot(1, 1, 1) 
    panel.tick_params(axis='both', which='major', labelsize=12)
    pyplot.plot(components**2, col)
    ##
    panel.xaxis.set_major_locator(majorXLocator)
    panel.xaxis.set_major_formatter(majorXFormatter)
    panel.yaxis.set_major_locator(majorYLocator)
    ##
    ##
    ## pyplot.annotate(bas_state, xy=(300, 0.07))
    pyplot.xlabel("Local Basis States")
    pyplot.ylabel("Squared Eigenvector Components")
    panel.xaxis.label.set_fontsize(18)
    panel.yaxis.label.set_fontsize(18)

