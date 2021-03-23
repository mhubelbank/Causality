import sys
import rv
import cGraph
from synth import getData
import independence
from Probability import Prob
from synth import synthDataGen
from standardize import standardize
import time

METHOD = 'prob'
#METHOD = 'fcit'
POWER = 1

args = sys.argv
test = 'Probability/Test/models/indCalibrationDat.csv'

r = getData.DataReader(test)
dat = r.read()
vars = dat.keys()
for var in vars:
    dat[var] = standardize(dat[var])
#print('dat = ', dat)
ps = Prob.ProbSpace(dat, power = 1, density = 1)

# List a variety of independent relationships
indeps = [('L1', 'L2'),
            ('L2', 'L3'),
            ('L1', 'L3'),
            ('E1', 'E2'),
            ('N1', 'N2'),
            ('L4', 'L5'),
            ('L5', 'L6'),
            ('L4', 'N3'),
            ('B', 'D', ['A']),
            ('A', 'C', ['B', 'D']),
            ('C', 'E2'),
            ('L6', 'L7', ['L3']),
            ('L4', 'L6', ['L3']),
            ('M1', 'E2'),
            ('M1', 'E2'),
            ]

# List a varieety of dependent relationships
deps = [('L3', 'L4'),
        ('L5', 'L2'),
        ('L6', 'L3'),
        ('L6', 'L7'),
        ('L7', 'L4'),
        ('E3', 'E1'),
        ('E3', 'E2'),
        ('M1', 'N2'),
        ('B', 'D', ['A', 'C']),
        ('B', 'A', 'C'),
        ('B', 'A', ['C', 'D']),
        ('C', 'B', 'A'),
        ('N1', 'N2', ['N3']),
        ('N3', 'E1', ['M1']),
        ]
print('Testing: ', test)
start = time.time()
minDep = ps.dependence('B', 'A', ['C', 'D'])
maxDep = ps.dependence('L5', 'L2')
maxIndep = ps.dependence('A', 'C', ['B', 'D'])
# Dependence Threshold Analysis
print('minDep = ', minDep)
print('maxDep = ', maxDep)
print('maxIndp = ', maxIndep)
print('Best Low Threshold = ', (minDep + maxIndep)/2)
print('Best High Threshold = ', maxDep + .01)
testVal = 0
condTestVal = 0
delta = .1

minIndep = 999999.0
maxDep = 0.0
cumIndep = 0.0
cumDep = 0.0
print()
print('Testing expected independents:')
for ind in indeps:
    if len(ind) == 2:
        x, y = ind
        z = []
    elif len(ind) == 3:
        x, y, z = ind
    else:
        print('*** Error, improperly specified independence =', ind)
    #xD = [dat[x]]
    #yD = [dat[y]]
    #zD = [dat[zvar] for zvar in z]
    pval = independence.test(ps, [x], [y], z, method = METHOD, power = POWER)
    #pval = ps.independence(x, y, z)
    if pval < minIndep:
        minIndep = pval
    print('independence', ind, '= ', pval)

print()
print('Testing expected dependents:')
print()
for dep in deps:
    if len(dep) == 2:
        x, y = dep
        z = []
    elif len(dep) == 3:
        x, y, z = dep
    else:
        print('*** Error, improperly specified independence =', dep)
    #xD = [dat[x]]
    #yD = [dat[y]]
    #zD = [dat[zvar] for zvar in z]
    pval = independence.test(ps, [x], [y], z, method = METHOD, power = POWER)
    #pval = ps.independence(x, y, z)
    if pval > maxDep:
        maxDep = pval
    print('independence', dep, ' = ', pval)
print()
print('Minimum value for expected independents = ', minIndep)
print('Maximum value for expected dependents =', maxDep)
print('Margin = ', minIndep - maxDep, '.  Positive margin is good.')
print('best threshold is: ', (minIndep + maxDep) / 2.0)
print()
end = time.time()
duration = end - start
print('Test Time = ', round(duration))