def hist_components_LDOS(eigvec_file, basis_state, energy_file, ExcE = True, N_value = 1, bins = 20, OutputCentroid = False, fig=False):
    ###
    '''LDOS Plot: Histogram of the squared components of a basis state as a function of
       the system eigenvectors (min basis state value = 0) and as a function of the system eigenvalues.

       Options: 
           ExcE    : Transform energy_file to excitation energies. Default value True
           N_value : renormalize by the system size (N_value = N)
           bins    : number of bins in the histogram
           OutputCentroid : If True outputs the centroid value for the bardiagram.
           fig     : plot bar fig

           '''
    #
    import numpy as np
    from matplotlib import pyplot
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    #
    from survival_probability_basis import read_eigenstates, read_energy_data
    #
    #
    eigstates=read_eigenstates(eigvec_file)
    ##
    eigenvalues = read_energy_data(energy_file)
    if (ExcE):
        E0 = eigenvalues[0]
        eigenvalues = eigenvalues - E0
    #
    eigenvalues = eigenvalues/N_value
    ## Energy bins
    energy_bins=np.linspace(eigenvalues[0], eigenvalues[-1],num = bins)
    col_values = np.zeros(len(energy_bins)) # Initialize histogram values
    #
    ## Basis state components
    components = eigstates[basis_state,:]**2   
    #
    columns = np.zeros(bins)
    ###
    idx = 0
    for index_bin in range(1,bins):
        #
        for aval_index in range(idx, len(eigenvalues)+1):
            #
            energy = eigenvalues[aval_index]
            #
            if (energy < energy_bins[index_bin]):
                columns[index_bin] = columns[index_bin] + components[aval_index]
                idx = idx + 1
            else:
                break
    #
    # Define bar centroids
    wdth = (energy_bins[1]-energy_bins[0])
    lngth = bins - 1
    centroids = energy_bins + wdth/2
    ##
    if (fig):
        figure = pyplot.figure()
        ax = pyplot.subplot(111)
        ax.bar(centroids[:lngth], columns[:lngth],width=wdth,alpha=0.5)
    ##
    if (OutputCentroid):
        return centroids[:lngth], columns[:lngth]
    else:
        output_data = np.zeros([2*bins+2,2])
        output_data[0,0] = energy_bins[0]
        output_data[0,1] = columns[0]
        #
        icount = 1
        for index_bin in range(1,bins):
            output_data[icount,0] = energy_bins[index_bin]
            output_data[icount,1] = columns[index_bin-1]
            output_data[icount+1,0] = energy_bins[index_bin]
            output_data[icount+1,1] = columns[index_bin]
            icount = icount + 2
        #   
        output_data[icount+1,0] = energy_bins[-1]
        output_data[icount+1,1] = columns[-1]
        return output_data
