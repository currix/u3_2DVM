def residuals_U3(U3H_FIT_param, input_filename, bin_path, BENT, exp_data_file, N_val, LMAX, VMAX, EMINL):
    import numpy as np
    from subprocess import Popen, PIPE
    #
    P11 = U3H_FIT_param['P11'].value
    P21 = U3H_FIT_param['P21'].value
    P22 = U3H_FIT_param['P22'].value
    P23 = U3H_FIT_param['P23'].value
    P31 = U3H_FIT_param['P31'].value
    P32 = U3H_FIT_param['P32'].value
    P33 = U3H_FIT_param['P33'].value
    P41 = U3H_FIT_param['P41'].value
    P42 = U3H_FIT_param['P42'].value
    P43 = U3H_FIT_param['P43'].value
    P44 = U3H_FIT_param['P44'].value
    P45 = U3H_FIT_param['P45'].value
    P46 = U3H_FIT_param['P46'].value
    P47 = U3H_FIT_param['P47'].value
    #
    # Prepare INPUT file
    input_f = open(input_filename, "w")
    output = '&INP0 BENT={}, exp_data_file="{}"  /\n'.format(BENT,exp_data_file)
    input_f.write(output)
    output = '&INP1 N_val={}, LMAX={}, VMAX={}, EMINL=.F. /\n'.format(N_val, LMAX, VMAX)
    input_f.write(output)
    output = '&INP2 IPRINT=0, DIS_RES = .T. /\n'
    input_f.write(output)
    output = '&INP1b P11={} /\n'.format(P11)
    # print output                # 
    input_f.write(output)
    output = '&INP2b P21={}, P22={}, P23={} /\n'.format(P21, P22, P23)
    # print output                # 
    input_f.write(output)
    output = '&INP3b P31={}, P32={}, P33={} /\n'.format(P31, P32, P33)
    input_f.write(output)
    output = '&INP4b P41={}, P42={}, P43={}, P44={}, P45={}, P46={}, P47={} /\n'.format(P41, P42, P43, P44, P45, P46, P47)
    input_f.write(output)
    input_f.close()
    #
    command = "echo " + input_filename + " | " + bin_path + "/chi2_U3_U2_gfortran"
    p = Popen([command], stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
    output, err = p.communicate()
    residuals = np.fromstring(output, dtype=float, sep='\n')
    #print "RESIDUALS ", residuals
    print np.sum(residuals**2)
    return residuals
