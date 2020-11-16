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
g = cGraph.cGraph(gnodes)

g.printGraph()

deps = g.computeDependencies(3)

#print('deps = ', deps)
g.printDependencies(deps)

# For dat file, use the input file name with the .csv extension
tokens = test.split('.')
testFileRoot = str.join('.',tokens[:-1])
datFileName = testFileRoot + '.csv'

d = getData.DataReader(datFileName)
data = d.read()
print('Testing Model for', len(deps), 'Independencies')
errors = 0
for dep in deps:
    node, isDep, x, y, zlist = dep
    X = data[x]
    Y = data[y]
    #print('X = ', X)
    #print('Y = ', Y)
    if zlist:
        Z = []
        for z in zlist:
            zdat = data[z]
            Z.append(zdat)
        fval = independence.test([X], [Y], Z)
    else:
        fval = independence.test([X], [Y])
    if isDep and fval > .1:
        print('*** Model Error(1): ', dep, fval)
        errors += 1
    elif not isDep and fval < .1:
        print('*** Model Error(2): ', dep, fval)
        errors += 1
    print('.',)
print('Model Testing Completed with', errors, 'errors.')


#print('data = ', data)