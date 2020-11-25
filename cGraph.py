import networkx
import math
import independence
VERBOSE = 1
# rvList is a list of random variable (rv) objects
# data is a dictionary keyed by random variable name, containing a list of observed values.
#   Each variable's list should be the same length
class cGraph:
    def __init__(self, rvList, data={}):
        self.g = networkx.DiGraph()
        self.rvDict = {}
        for rv in rvList:
            if rv.name in self.rvDict.keys():
                raise 'Duplicate variable name = ' + rv.name
            self.rvDict[rv.name] = rv
        self.g.add_nodes_from(self.rvDict.keys())
        edges = []
        for rv in rvList:
            for pa in rv.parentNames:
                edges.append((pa, rv.name))
        self.g.add_edges_from(edges)
        self.data = data
        edges = self.g.edges()
        self.edgeDict = {}
        for edge in edges:
            s, d = edge
            if s in self.edgeDict.keys():
                self.edgeDict[s].append(edge)
            else:
                self.edgeDict[s] = [edge]
            if d in self.edgeDict.keys():
                self.edgeDict[d].append(edge)
            else:
                self.edgeDict[d] = [edge]
        
    def printGraph(self):
        print('Nodes:', self.g.nodes())
        print('Edges:', self.g.edges())


    def getAdjacencies(self, node):
        return self.edgeDict[node]

    def combinations(self, inSet):
        c = []
        
        for i in range(len(inSet)):
            u = inSet[i]
            for v in inSet[i+1:]:
                c.append((u, v))
        return c

    def makeDependency(self, u, v, w, isDep):
        if u < v:
            d = (u, v, w, isDep)
        else:
            d = (v, u, w, isDep)
        return d

    def calcNDependencies(self, order, n=0):
        def combinations(n, r):
            return math.factorial(n)/ (math.factorial(r) * math.factorial(n - r))
        if n == 0:
            n = len(g.nodes())
        nDeps = n * (n-1) / 2
        nCDeps = 0
        for o in range(order):
            r = o + 1
            if r > n - 2:
                break
            nCDeps += nDeps * combinations(n-2, r)
        return (nDeps, nCDeps, nDeps + nCDeps)

    def getCombinations(self, nodes, order):
        from itertools import combinations
        nodes = self.g.nodes()
        allCombos = []
        for o in range(1, order):
            combos = combinations(nodes, o)
            allCombos += combos
        return allCombos


    def computeDependencies(self, order):
        deps = []
        nodes = list(self.g.nodes())
        print('nodes = ', nodes)
        cNodes = self.getCombinations(nodes, order)
        for i in range(len(nodes)):
            node1 = nodes[i]
            if not self.rvDict[node1].isObserved:
                continue
            for j in range(i, len(nodes)):
                node2 = nodes[j]
                if node1 == node2 or not self.rvDict[node2].isObserved:
                    continue
                isSeparated = networkx.d_separated(self.g, {node1}, {node2}, {})
                dep = self.makeDependency(node1, node2, None, not isSeparated)
                deps.append(dep)
                for c in cNodes:
                    if node1 in c or node2 in c:
                        continue
                    # Verify that every member of c is observed.  If not, we skip this combo.
                    allObserved = True
                    for m in c:
                        if not self.rvDict[m].isObserved:
                            allObserved = False
                            break
                    if not allObserved:
                        continue
                    isSeparated = networkx.d_separated(self.g, {node1}, {node2}, set(c))
                    dep = self.makeDependency(node1, node2, c, not isSeparated)
                    deps.append(dep)
        return deps
    
    def formatDependency(self, dep):
        # dep is:  from, to, given, isDependent
        u, v, w, isDep = dep
        if isDep:
            rel = 'is not independent from'
        else:
            rel = 'is independent from'
        if w is None:
            given = ''
        else:
            given = 'given ' + str(w)
        out = u + ' ' + rel + ' ' + v + ' ' + given
        return out

    def printDependencies(self, deps):
        print('Implied Dependencies:\n')
        for d in deps:
            print(self.formatDependency(d))

    # Test the model for consistency with a set of data.
    # Format for data is {variableName: [variable value]}.
    # That is, a dictionary keyed by variable name, containing
    # a list of data values for that variable's series.
    # The lengths of all variable's lists should match.
    # That is, the number of samples for each variable must
    # be the same. 
    # Returns [(errType, x, y, z, isDep, errStr)]
    # Where errType = 1 (Expected ndependence not observed) or 
    #                   2 (Expected dependence not observed)
    #       x, y, z are each a list of variable names that
    #               comprise the statement x _||_ y | z.
    #               That is x is independent of y given z.
    #       isDep True if a dependence is expected.  False for 
    #               independence
    #       pval -- The p-val returned from the independence test
    #       errStr A human readable error string describing the error
    #
    def TestModel(self, data, order=3):
        errors = []
        deps = self.computeDependencies(order)
        if VERBOSE:
            print('Testing Model for', len(deps), 'Independencies')
        for dep in deps:
            x, y, zlist, isDep = dep
            X = data[x]
            Y = data[y]
            #print('X = ', X)
            #print('Y = ', Y)
            Z = []
            if zlist:
                for z in zlist:
                    zdat = data[z]
                    Z.append(zdat)
                pval = independence.test([X], [Y], Z)
            else:
                pval = independence.test([X], [Y])
            errType = 0
            if isDep and pval > .1:
                errStr = 'Warning (Type 2) -- Expected: ' +  self.formatDependency(dep) + ' but no dependence detected.  P-val = ' + str(pval)
                errType = 2
            elif not isDep and pval < .1:
                errStr = 'Error (Type 1) -- Expected: ' + self.formatDependency(dep) + ' but dependence was detected. P-val = ' + str(pval)
                errType = 1
            if errType:
                if VERBOSE:
                    print('***', errStr)
                errors.append((errType, [x], [y], list(zlist), isDep, pval, errStr))
            elif VERBOSE:
                print('.',)
        if VERBOSE:
            print('Model Testing Completed with', len(errors), 'error(s).')
        return errors
        


    