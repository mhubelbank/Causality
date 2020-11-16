import networkx
import math

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

    def makeDependency(self, node, isDep, u, v, w):
        if u < v:
            d = (node, isDep, u, v, w)
        else:
            d = (node, isDep, v, u, w)
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
            for j in range(i, len(nodes)):
                node2 = nodes[j]
                if node1 == node2:
                    continue
                isSeparated = networkx.d_separated(self.g, {node1}, {node2}, {})
                dep = self.makeDependency(None, not isSeparated, node1, node2, None)
                deps.append(dep)
                for c in cNodes:
                    if node1 in c or node2 in c:
                        continue
                    isSeparated = networkx.d_separated(self.g, {node1}, {node2}, set(c))
                    dep = self.makeDependency(None, not isSeparated, node1, node2, c)
                    deps.append(dep)
        return deps
                   
                    
    def computeDependenciesOld(self):
        # dependencies := [(isDependent, u, v, w)]
        deps = []
        for node1 in self.g.nodes():
            edges = self.getAdjacencies(node)
            inbound = []
            outbound = []
            for edge in edges:
                s, d = edge
                if s == node:
                    outbound.append(d)
                else:
                    inbound.append(s)
            if inbound and not outbound:
                combos = self.combinations(inbound)
                for neighbor in inbound:
                    deps.append(self.makeDependency(node, True, node, neighbor, None))
                for c in combos:
                    u, v = c
                    deps.append(self.makeDependency(node, False, u, v, None))
                    deps.append(self.makeDependency(node, True, u, v, node))
            elif outbound and not inbound:
                combos = self.combinations(outbound)
                for neighbor in outbound:
                    deps.append(self.makeDependency(node, True, node, neighbor, None))
                for c in combos:
                    u, v = c
                    deps.append(self.makeDependency(node, True, u, v, None))
                    deps.append(self.makeDependency(node, False, u, v, node))
            elif outbound and inbound:
                for neighbor in inbound:
                    deps.append(self.makeDependency(node, True, node, neighbor, None))
                for neighbor in outbound:
                    deps.append(self.makeDependency(node, True, node, neighbor, None))
                iCombos = self.combinations(inbound)
                for c in iCombos:
                    u, v = c
                    deps.append(self.makeDependency(node, False, u, v, None))
                    deps.append(self.makeDependency(node, True, u, v, node))
                oCombos = self.combinations(outbound)
                for c in oCombos:
                    u, v = c
                    deps.append(self.makeDependency(node, True, u, v, None))
                    deps.append(self.makeDependency(node, False, u, v, node))
                for i in inbound:
                    for j in outbound:
                        deps.append(self.makeDependency(node, False, i, j, node))
            else:
                # No neighbors
                pass
            # Now reduce the local dependency lists
            depDict = {}
            condDepDict = {}
            indDict = {}
            condIndDict = {}
            outDeps = []
            for d in deps:
                node, isDep, u, v, w = d
                if w is None:
                    # Non conditional
                    if isDep:
                        depDict[(u, v)] = d
                    else:
                        indDict[(u, v)] = d
                else:
                    if isDep:
                        condDepDict[(u, v, w)] = d
                    else:
                        condIndDict[(u, v, w)] = d
            # Dedup any recods by reconstituting from dictoinary
            deps = list(depDict.values()) + list(indDict.values()) + list(condDepDict.values()) + list(condIndDict.values())
            for d in deps:
                node, isDep, u, v, w = d
                if w is None:
                    # Non conditional
                    if isDep:
                        outDeps.append(d)
                    elif (u, v) not in depDict.keys():
                        outDeps.append(d)
                else:
                    if isDep:
                        outDeps.append(d)
                        #if (u, v) in depDict.keys():
                            # If there is a direct dependency between u and v, then invert the conditional rule
                            #d = (node, False, u, v, w)
                    elif (u, v, w) not in condDepDict.keys():
                        #if (u, v) in depDict.keys():
                            # If there is a direct dependency between u and v, then invert the conditional rule
                        #    d = (node, True, u, v, w)
                        outDeps.append(d)
        
        return outDeps

    def printDependencies(self, deps):
        print('Implied Dependencies:\n')
        for d in deps:
            node, isDep, u, v, w = d
            if isDep:
                rel = 'is not independent from'
            else:
                rel = 'is independent from'
            if w is None:
                given = ''
            else:
                given = 'given ' + str(w)
            
            print(u, rel, v, given)

    def testModel(self):
        d = self.computeDependencies()
        return d      
        


    