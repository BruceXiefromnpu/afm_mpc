# call_pc.py
def place_varg(A, B, Ns, Nu, P_real, P_imag):
    # import numpy.distutils.system_info as sysinfo
    # p = sysinfo.get_info('lapack_opt')
    # print(p)

    import numpy as np
    #np.show_config()
    # import numpy as np
    # import scipy.linalg as linalg
    import sys
    # #import ipdb
    # # old_path = '/usr/local/lib/python2.7/dist-packages/control-dev-py2.7.egg'
    # # paths = sys.path
    # # if old_path in paths:
    # #     sys.path.remove(old_path)
    slycot_path = '/home/arnold/pythonBox/control_dev/slycot-rabraker/build/lib.linux-x86_64-2.7/'
    if slycot_path not in sys.path:
        sys.path.append(slycot_path)

    # # A = np.array([[0, 1],[100, 0]])
    # # B = np.array([[0],[1]])

    # # Pdes = np.array([-20 + 10*1j, -20 - 10*1j])

    # # #######################################################
    from slycot import sb01bd

    A = np.array(A)
    A = A.reshape((Ns, Ns), order='F')
    # A = np.array([[0.9904, -0.0897,    0.0399],
    #               [0.1594,    0.9928,    0.0032],
    #               [0,         0,    1.0000]])
    B = np.array(B, order='F')
    B = B.reshape((Ns, Nu))
    Pdes = np.array(P_real) + np.array(P_imag)*1j

    # # # SB01BD sets eigenvalues with real part less than alpha
    # # # We want to place all poles of the system => set alpha to minimum
    #alpha = min(np.linalg.eigvals(A).real)*0
    alpha = 0.1
    # # # # Call SLICOT routine to place the eigenvalues
    A_z,w,nfp,nap,nup,K,Z = sb01bd(B.shape[0], B.shape[1],
                               len(Pdes), alpha, A, B, Pdes, 'D')
    #print 'Placed Poles: %s'%np.linalg.eigvals(A + B.dot(K))

    import array
    return array.array('d', K.flatten())


#import array
# p_imag = array.array('d', [0.0, 0.05704602973441748, -0.05704602973441748])

# p_real = array.array('d', [0.8, 0.9012258993608832, 0.9012258993608832])

# b = array.array('d', [0.03985622250751898, 0.0031936070602840744, 1.0])
# a = array.array('d', [0.9904230107639097, 0.15942489003007596, 0.0,
#                 -0.0896765006419177, 0.9928143841143608, 0.0,
#                 0.03985622250751898, 0.0031936070602840744, 1.0])

# K = place_varg(a, b, 3, 1, p_real, p_imag)

# print(K)
