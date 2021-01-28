from Statistics import Stats
from synth import getData
import sys

def run(filename):
    r = getData.DataReader(filename)
    dat = r.read()
    stats = Stats.Stats(dat)
    print ('discretized = ', stats.discretized)
    print ('discSpec = ', stats.discSpecs)
    print('pdf(A) = ', stats.pdf('A'))
    print('pdf(C) = ', stats.pdf('C'))
    print('pdf(B) = ', stats.pdf('B'))

    print('P(B=0|A=1) = ', stats.condProb('B', 0, ('A', 1)))
    print('P(B=1|A=1) = ', stats.condProb('B', 1, ('A', 1)))
    print('P(B=2|A=1) = ', stats.condProb('B', 2, ('A', 1)))
    print('P(B=0|A=0) = ', stats.condProb('B', 0, ('A', 0)))
    print('P(B=1|A=0) = ', stats.condProb('B', 1, ('A', 0)))
    print('P(B=2|A=0) = ', stats.condProb('B', 2, ('A', 0)))
    print('P(B=0) = ', stats.prob('B', 0))
    print('P(B=1) = ', stats.prob('B', 1))
    print('P(B=2) = ', stats.prob('B', 2))
  
if __name__ == '__main__':
    filename = None
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        run(filename)
    else:
        print ('first parameter is the SEM filename -- Required')
