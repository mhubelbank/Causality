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
    print('stats(A) =  ', samp.fieldStats('A'))
    print('stats(C) = ', samp.fieldStats('C'))
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
    print('E(C|B=0, A=1) =', samp.distr('C', [('B', 0), ('A', 1)]).E())
    print('E(C|B=0, A=2) =', samp.distr('C', [('B', 0), ('A', 2)]).E())
    print('E(C|B=1, A=1) =', samp.distr('C', [('B', 1), ('A', 1)]).E())
    print('E(C|B=1, A=2) =', samp.distr('C', [('B', 1), ('A', 2)]).E())
    print('E(C|A=1) =', samp.distr('C', [('A', 1)]).E())
    print('E(C|A=2) =', samp.distr('C', [('A', 2)]).E())
    print('E(C|A=1,B) = ', samp.distr('C' , [('A', 1), 'B']).E())
    print('E(C|A=2,B) = ', samp.distr('C' , [('A', 2), 'B']).E())
    print('ACE(A,C) = ', samp.distr('C' , [('A', 2), 'B']).E() - samp.distr('C' , [('A', 1), 'B']).E())

if __name__ == '__main__':
    filename = None
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        run(filename)
    else:
        print ('first parameter is the SEM filename -- Required')
