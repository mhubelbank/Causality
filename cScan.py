import sys
import rv
import cGraph
from synth import getData
import independence
import time

args = sys.argv
if (len(args) > 1):
    test = args[1]

f = open(test, 'r')
exec(f.read(), globals())

print('Testing: ', test, '--', testDescript)
print()
start = time.time()
gnodes = []

# 'model' is set when the text file is exec'ed
for var in model:
    observed = True
    dType = 'Numeric'
    if type(var) == type((0,)):
        name, parents = var[:2]
        if len(var) >= 3:
            observed = var[2]
        if len(var) >= 4:
            dType = var[3]
    else:
        name = var
        parents = []
    gnode = rv.rv(name, parents, observed, dType, None, None)
    gnodes.append(gnode)

# For dat file, use the input file name with the .csv extension
tokens = test.split('.')
testFileRoot = str.join('.',tokens[:-1])
datFileName = testFileRoot + '.csv'

d = getData.DataReader(datFileName)
data = d.read()

g = cGraph.cGraph(gnodes, data)

#g.printGraph()

exos = g.findExogenous()
print()
print('Exogenous variables = ', exos)
cvMap = g.findChildVars(exos)
print('CVMap = \n', cvMap)
print()
dm = {} # Discovered Model
for t in cvMap:
    parent, child = t
    if child not in dm.keys():
        dm[child] = []
    if parent not in dm.keys():
        dm[parent] = []
    dm[child].append(parent)
dm2 = []
for var in dm.keys():
    parents = dm[var]
    parents.sort()
    dm2.append((len(parents), var, parents))
dm2.sort()
dm3 = [varSpec[1:] for varSpec in dm2]

print('Discovered Model = \n', dm3, '\n')
print()

clusters = {}
for var in exos:
    clusters[(var,)] = [var]
for spec in dm3:
    var, parents = spec
    if var in exos:
        continue
    clustKey = tuple(parents)
    if clustKey not in clusters.keys():
        clusters[clustKey] = []
    clusters[clustKey].append(var)
print('clusters = ', clusters)
print()
clusterLinks = []
for clustId in clusters.keys():
    clustData = {}
    clustGNodes = []
    clustVars = clusters[clustId]
    memberCnt = 0
    for var in clustVars:
        clustData[var] = data[var]
        clustGNodes.append(rv.rv(var, [], True, dType, None, None))
        if var not in exos:
            memberCnt += 1
    #if memberCnt > 1:
    if len(clustVars) > 1:
        #print('clustGNodes = ', [rv.name for rv in clustGNodes])
        cg = cGraph.cGraph(clustGNodes, clustData)
        cOrder = cg.causalOrder()
        print('order for cluster', clustId, ' = ', cOrder)
        clusters[clustId] = cOrder
        for i  in  range(1, len(cOrder)):
            testVar = cOrder[i]
            parent = cOrder[i-1]
            for exo in list(clustId):
                exDep1 = g.iProb.dependence(parent, exo)
                exDep2 = g.iProb.dependence(testVar, exo)
                epsilon = .05
                if exDep2 - exDep1 > epsilon:
                    print(testVar, 'contains more', exo, 'than ', parent, exDep1, exDep2, '.  This variable must be independently connected to Cluster = ', (exo,))
                    clusterLinks.append(((exo,), testVar))
    cOrder = clusters[clustId]
    if len(clustId) > 1:
        for exo in list(clustId):
            clusterLinks.append(((exo,), cOrder[0]))
    else:
        clustExo = clustId[0]
        for var in cOrder:
            if var == clustExo:
                continue
            clusterLinks.append((exo, var))
    
    clusterLinks2 = clusterLinks
    clusterLinks = []
    for link in clusterLinks2:
        parent, child = link
        if type(parent) == type((1,)):
            # Parent is a cluster.  Try to resolve it.
            cOrder = clusters[parent]
            if len(cOrder) == 1:
                clusterLinks.append((cOrder[0], child))
            else:
                # Add more resolution here.
                clusterLinks.append(link)
        else:
            # Resolved.
            clusterLinks.append(link)
        

            

print('clusterLinks = ', clusterLinks)

                    





end = time.time()
duration = end - start
print('Test Time = ', round(duration))