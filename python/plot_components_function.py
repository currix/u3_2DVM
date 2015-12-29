def plot_components(eigvec_file, eigstate_index, save_data = False, X_vector = 0, file_name = "sqrd_component_out.dat", majorX = 50, majorY = 0.02, col="b-o"):
    ###
    '''Plot or save the squared components of the eigstate_index-th eigenvector (min value = 0).

       Options: 
           save_data : if True do no plot the data but save them.
           X_vector  : plot as a function of the X_vector instead of simple indices.
           file_name  : squared component filename if save_data = True 

           '''
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
    #
    majorYLocator   = MultipleLocator(majorY)
    #
    components = eigstates[:,eigstate_index]**2    ## Eigenstate's squared components
    ##
    if (save_data):
        # Saving data
        if (isinstance(X_vector,np.ndarray)):
            data_save = np.transpose(np.array([X_vector,components])) 
        else:
            data_save = np.transpose(np.array([range(len(components)),components])) 
            ##
        np.savetxt(file_name, data_save)
    else:
        # Plot data
        panel = pyplot.subplot(1, 1, 1) 
    ##
    ##
        if (not isinstance(X_vector,np.ndarray)):
            majorXLocator   = MultipleLocator(majorX)
            majorXFormatter = FormatStrFormatter('%d')
            panel.tick_params(axis='both', which='major', labelsize=12)
            pyplot.plot(components, col)
            panel.xaxis.set_major_locator(majorXLocator)
            panel.xaxis.set_major_formatter(majorXFormatter)
        else:
            panel.tick_params(axis='both', which='major', labelsize=12)
            pyplot.plot(X_vector, components, col, markersize = 6)        
    ##
            panel.yaxis.set_major_locator(majorYLocator)
    ##
    ##
            panel.text(0.7, 0.03, r'$k = 49$', fontsize=15)
            pyplot.plot([0.4166666666666667, 0.4166666666666667], [0, 0.039], 'k-.', lw=2)
    ## pyplot.annotate("$k = 148$", xy=(0.1, 0.025))
    ## pyplot.xlabel("Normalized Basis State Energy $e_\omega$")  ## SO(4) basis
            pyplot.xlabel("Normalized Basis State Energy $e_n$") ## U(3) basis
    ## pyplot.ylabel("$|C_\omega^{(k)}|^2$") ## SO(4) basis
            pyplot.ylabel("$|C_n^{(k)}|^2$") ## U(3) basis
            panel.xaxis.label.set_fontsize(18)
            panel.yaxis.label.set_fontsize(18)
    return
            
