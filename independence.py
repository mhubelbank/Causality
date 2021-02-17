from Probability import Prob

# Tests for independence between two sets of variables,
# optionally given a third set.
# X, Y, and Z are each a list of data series', one series
# per variable, with each series representing the values at
# each sample.
# For example, Z with 2 variables and N samples would be:
# [[v1[0], ... v1[N-1]], [v2[1], ... , v2[N-1]]]
# Returns a p-value for the null hypothesis that:
# X and Y are Dependent given Z.  P-values less than
# .05 provide a 95% confidence that the variables are dependent.
# P-values > .05 imply independence (i.e. lack of proof of dependence).
def testFCIT(X, Y, Z=[]):
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
        pval = fcit.test(Xa, Ya, num_perm = 100, prop_test = .40)
    return pval

def testSDCIT(X, Y, Z=[]):
    import numpy as np
    from sdcit.sdcit_mod import SDCIT
    from sdcit.utils import rbf_kernel_median

    Xa = np.array(X).transpose()
    Ya = np.array(Y).transpose()
    if not Z:
        return testFCIT(X, Y)
    #Za = np.zeros((len(Z[0]),))
    #for i in range(len(Z)):
    #    Za = np.array(Z[i]) + Za
    Za = np.array(Z).transpose()
    Kx, Ky, Kz = rbf_kernel_median(Xa, Ya, Za)
    test_stat, p_value = SDCIT(Kx, Ky, Kz)
    #print('p = ', p_value)
    return p_value

def testProb(X, Y, Z=[], level=2):
    X = X[0]
    Y = Y[0]
    d = {}
    d['X'] = X
    d['Y'] = Y
    zNames = []
    for i in range(len(Z)):
        varName = 'Z' + str(i)
        d[varName] = Z[i]
        zNames.append(varName)
    s = Prob.Sample(d, .2)
    dep = s.dependence('X', 'Y', zNames, level=level)
    print('dep = ', dep)
    if dep <= .15:
        return 1.0
    else:
        return 0

def test(X, Y, Z=[]):
    #p_val = testFCIT(X, Y, Z)
    p_val = testProb(X, Y, Z, level = 2)
    #p_val = testSDCIT(X, Y, Z)
    return p_val