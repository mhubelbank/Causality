import sys
import rv
import cGraph
from synth import getData
import independence

args = sys.argv
if (len(args) > 1):
    test = args[1]

f = open(test, 'r')
exec(f.read(), globals())

print('Testing: ', test, '--', testDescript)
gnodes = []
# 'model' is set when the text file is exec'ed
for var in model:
    observed = True
    dType = 'Numeric'
    name, parents = var[:2]
    if len(var) >= 3:
        observed = var[2]
    if len(var) >= 4:
        dType = var[3]
    gnode = rv.rv(name, parents, observed, dType, None, None)
    gnodes.append(gnode)

# For dat file, use the input file name with the .csv extension
tokens = test.split('.')
testFileRoot = str.join('.',tokens[:-1])
datFileName = testFileRoot + '.csv'

d = getData.DataReader(datFileName)
data = d.read()

g = cGraph.cGraph(gnodes, data)

g.printGraph()

deps = g.computeDependencies(3)

#print('deps = ', deps)
g.printDependencies(deps)
print('P(C)')
print('E(C | do(A=1)) = ', g.intervene('C', [('A', 1)]).E())
print('E(C | do(A=2)) = ', g.intervene('C', [('A', 2)]).E())
print('E(C | A=1, B) = ', g.prob.distr('C', [('A', 1), 'B']).E())
print('E(C | A=2, B) = ', g.prob.distr('C', [('A', 2), 'B']).E())
print('ACE(A, C) = ', g.ACE('A', 'C'))
print('ACE(C, A) = ', g.ACE('C', 'A'))

#results = g.TestModel()
#conf = results[0]
#print('\nConfidence = ', round(conf * 100, 1), '%, testPerType = ', results[2], ', errsPerType = ', results[3])
#print('TestResults = ', results)
