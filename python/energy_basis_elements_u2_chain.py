def energy_basis_u2(N_value, n_value, l_value, xi_value, alpha_value = 0.0):
    #
    '''Computes the diagonal matrix element <[N] n L|H(xi,alpha)|[N] n L> '''
    #
    return (1-xi_value) * n_value + (alpha_value/(N_value - 1))*n_value*(n_value+1.0) + (xi_value/(N_value-1)) * (N_value*(N_value + 1) - (N_value - n_value)*(n_value + 2.0) - (N_value - n_value + 1.0)*n_value - l_value*l_value)
