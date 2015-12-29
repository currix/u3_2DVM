def energy_basis_so3(N_value, omega_value, l_value, xi_value, alpha_value = 0.0):
    #
    '''Computes the diagonal matrix element <[N] w L|H(xi, alpha)|[N] w L> '''
    #
    import numpy as np
    #
    # U(2) number operator n contribution
    ndiagmatel = (N_value - omega_value)*((omega_value - l_value + 2.0)*(omega_value - l_value + 1.0) + (omega_value + l_value + 2.0)*(omega_value + l_value + 1.0))/(2.0*(2.0*omega_value+1.0)*(2.0*omega_value+3.0)) +  (N_value + omega_value + 1.0)*((omega_value+l_value)*(omega_value+l_value-1.0) + (omega_value-l_value)*(omega_value-l_value-1.0)) / (2.0*(2.0*omega_value-1.0)*(2.0*omega_value+1.0))
    ##
    n2diagmatel = 0.0
    if (alpha_value != 0.0):
        #
        n2diagmatel = ndiagmatel*(ndiagmatel + 1.0) + B0(N_value, omega_value, l_value)**2 + B0(N_value, omega_value-2, l_value)**2
    ##
    # SO(3) casimir operator contribution
    pairdiagmatel = N_value*(N_value + 1.0) - omega_value*(omega_value + 1.0)
    #
    # Note the N-1 normalization in the SO(3) casimir operator and the - sign
    return (1-xi_value)*ndiagmatel +  (alpha_value/(N_value - 1.0))*n2diagmatel + (xi_value/(N_value - 1.0))*pairdiagmatel
########
# Anharmonicity contribution n(n+1)
def B0(N_value, omega, l_value):
    import numpy as np
    #
    B0value = np.sqrt((N_value - omega)*(N_value + omega + 3.0)*(omega-l_value+2.0)*(omega+l_value+2.0)*(omega+l_value+1.0)*(omega-l_value+1.0)/((2.0*omega+1.0)*(2.0*omega+3.0)**2*(2.0*omega+5.0)) )
        #
    return B0value
#
