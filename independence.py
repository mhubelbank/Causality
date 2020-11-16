def test(X, Y, Z=[]):
    import numpy as np
    from fcit import fcit

    Xa = np.array(X).transpose()
    #print('xshape = ', Xa.shape)
    Ya = np.array(Y).transpose()
    if Z:
        Za = np.array(Z).transpose()
        #print('zshape = ', Za.shape)
        pval = fcit.test(Xa, Ya, Za)
    else:
        pval = fcit.test(Xa, Ya)
    return pval
