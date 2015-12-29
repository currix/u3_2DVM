def hist_components_LDOS(eigvec_file, basis_index, energy_file, N_value = 1, bins = 20, fig=False):
    ###
    '''LDOS Plot: Histogram of the squared components of the index-th basis state as a function of
       the system eigenvectors (min value = 0) and as a function of the system eigenvalues.

       Options: 

           bins    : number of bins in the histogram
           N_value : renormalize by the system size
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
    eigenvalues = read_energy_data(energy_file)/N_value
    #
    ## Energy bins
    energy_bins=np.linspace(eigenvalues[0], eigenvalues[-1],num = bins)
    col_values = np.zeros(len(energy_bins)) # Initialize histogram values
    #
    ## Basis state components
    components = eigstates[basis_index,:]**2   
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
    return centroids[:lngth], columns[:lngth]
