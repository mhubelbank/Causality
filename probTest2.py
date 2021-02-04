from Probability.Prob import Sample
from synth import getData
import sys

def run(filename):
    r = getData.DataReader(filename)
    dat = r.read()
    samp = Sample(dat, density=1)
    #print ('discretized = ', samp.discretized)
    #print ('discSpecs = ', samp.discSpecs)
    #print('jointVals = ', samp.jointCondSpecs(['B', 'A']))
    print('stats(B) =  ', samp.fieldStats('B'))
    print('stats(A) =  ', samp.fieldStats('A'))
    print('stats(C) = ', samp.fieldStats('C'))
    print('values(B) = ', samp.getMidpoints('B'))
    print('E(A) = ', samp.distr('A').E())
    print('E(C) = ', samp.distr('C').E())
    print('E(B) = ', samp.distr('B').E())
    print('P(B=0) = ', samp.P('B', 0))
    print('P(B=1) = ', samp.P('B', 1))
    print('P(B=2) = ', samp.P('B', 2))
    print('E(A) = ', samp.distr('A').E())
    print('E(A|B=0) =', samp.PD('A', ('B', 0)).E())
    print('E(A|B=1) =', samp.distr('A', ('B', 1)).E())
    print('E(C) = ', samp.distr('C').E())
    print('P(B=-1)', samp.prob('B', -1))
    print('P(B=0)', samp.prob('B', 0))
    print('P(B=0)', samp.prob('B', 1))

    print('P(A=1,B=2) = ', samp.jointProb([('A', 1), ('B', 2)]))
    print('P(A=2,B=1) = ', samp.jointProb([('A', 2), ('B', 1)]))
    print('P(B=1,A=2) = ', samp.jointProb([('B', 1), ('A', 2)]))

    print('E(C | B)', samp.distr('C', 'B').E())
    print('E(C)', samp.distr('C').E())
    bVals = samp.getMidpoints('B')
    aMean = samp.distr('A').E()
    aLow = aMean - .5
    aHigh = aMean + .5
    for bVal in bVals:
        print('E(C|A=', aLow,',B=', bVal, ') = ', samp.distr('C', [('A', aLow), ('B', bVal)]).E())
        print('E(C|A=', aHigh, ',B=', bVal, ') = ', samp.distr('C', [('A', aHigh), ('B', bVal)]).E())


    print('E(C|A=aLow) =', samp.distr('C', [('A', aLow)]).E())
    print('E(C|A=aHigh) =', samp.distr('C', [('A', aHigh)]).E())

    print('E(C|A=aLow,B) = ', samp.distr('C' , [('A', aLow), 'B']).E())
    print('E(C|A=aHigh,B) = ', samp.distr('C' , [('A', aHigh), 'B']).E())

    print('ACE(A,C) = ', samp.distr('C' , [('A', aHigh), 'B']).E() - samp.distr('C' , [('A', aLow), 'B']).E())

if __name__ == '__main__':
    filename = None
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        run(filename)
    else:
        print ('first parameter is the SEM filename -- Required')
