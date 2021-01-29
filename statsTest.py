from Statistics import Stats
from synth import getData
import sys

def run(filename):
    r = getData.DataReader(filename)
    dat = r.read()
    stats = Stats.Stats(dat)
    #print ('discretized = ', stats.discretized)
    #print ('discSpec = ', stats.discSpecs)
    print('E(A) = ', stats.expected(stats.pdf('A')))
    #print('P(C) = ', stats.pdf('C'))
    print('E(C) = ', stats.expected(stats.pdf('C')))
    print('E(B) = ', stats.expected(stats.pdf('B'), True))
    print('stats(A) =  ', stats.fieldStats('A'))
    print('stats(C) = ', stats.fieldStats('C'))
    #print('P(A|B=0) =', stats.condPDF('A', ('B', 0)))
    #print('P(A|B=1) =', stats.condPDF('A', ('B', 1)))
    print('E(A|B=0) =', stats.expected(stats.condPDF('A', ('B', 0))))
    print('E(A|B=1) =', stats.expected(stats.condPDF('A', ('B', 1))))
    #print('P(C|B=0, A=0) =', stats.expected(stats.condPDF('C', [('B', 0), ('A', 0)])))
    #print('P(C|B=0, A=0) =', stats.expected(stats.condPDF('C', [('B', 0), ('A', 1)])))
    print('E(C|B=0, A=1) =', stats.expected(stats.condPDF('C', [('B', 0), ('A', 1)])))
    print('E(C|B=0, A=2) =', stats.expected(stats.condPDF('C', [('B', 0), ('A', 2)])))
    print('E(C|B=1, A=1) =', stats.expected(stats.condPDF('C', [('B', 1), ('A', 1)])))
    print('E(C|B=1, A=2) =', stats.expected(stats.condPDF('C', [('B', 1), ('A', 2)])))
    print('E(C|A=1) =', stats.expected(stats.condPDF('C', [('A', 1)])))
    print('E(C|A=2) =', stats.expected(stats.condPDF('C', [('A', 2)])))
    print('E(C|A=1,B) = ', stats.expected(stats.condPDF('C' , [('A', 1), ('B', None)])))
    print('E(C|A=2,B) = ', stats.expected(stats.condPDF('C' , [('A', 2), ('B', None)])))
    print('ACE(A,C) = ', stats.expected(stats.condPDF('C' , [('A', 2), ('B', None)]))-stats.expected(stats.condPDF('C' , [('A', 1), ('B', None)])))
    #print('E(C|A=0,B) = ', stats.expected(stats.condPDF('C' , [('A', 0), ('B', None)]), True))
    #print('E(C|A=1,B) = ', stats.expected(stats.condPDF('C' , [('A', 1), ('B', None)]), True))
    # print('P(B=1|A=1) = ', stats.condProb('B', 1, ('A', 1)))
    # print('P(B=2|A=1) = ', stats.condProb('B', 2, ('A', 1)))
    # print('P(B=0|A=0) = ', stats.condProb('B', 0, ('A', 0)))
    # print('P(B=1|A=0) = ', stats.condProb('B', 1, ('A', 0)))
    # print('P(B=2|A=0) = ', stats.condProb('B', 2, ('A', 0)))
    # print('P(B=0) = ', stats.prob('B', 0))
    # print('P(B=1) = ', stats.prob('B', 1))
    # print('P(B=2) = ', stats.prob('B', 2))
    # print('P(B |A=0)', stats.condPDF('B', ('A', 1)))
    # print('P(A | B =1)', stats.condPDF('A', ('B', 1)))
    # print('P(B | B=1, C=1)', stats.condPDF('B', [('B', 1), ('C', 1)]))
    # print('P(C | B=0)', stats.condPDF('C', [('B', 0)]))
    # print('P(C | B=1)', stats.condPDF('C', [('B', 1)]))
    # print('P(C | B=2)', stats.condPDF('C', [('B', 2)]))
    # print('P(C | B=0, A)', stats.condPDF('C', [('B', 1), ('A', None)]))

if __name__ == '__main__':
    filename = None
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        run(filename)
    else:
        print ('first parameter is the SEM filename -- Required')
