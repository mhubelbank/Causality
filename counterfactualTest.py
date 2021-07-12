import sys
import rv
import counterfactual
from synth import getData

args = sys.argv
if (len(args) > 1):
    test = args[1]
    pdf = True if len(args) > 2 and args[1] else False

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
testFileRoot = str.join('.', tokens[:-1])
datFileName = testFileRoot + '.csv'

d = getData.DataReader(datFileName)
data = d.read()

cf = counterfactual.Counterfactual(gnodes, data, model, varEquations)  # SEM is set when the text file is exec'ed
g = cf.cGraph

sally = [('G', 0), ('H', 61), ('W', 110), ('S', 4), ('B', 22.78), ('T', 1), ('BR', 0.019), ('R', 0)]
print(cf.cf(sally, 'R', ('S', 3), pdf))
print(cf.cf(sally, 'R', ('S', 2), pdf))
print(cf.cf(sally, 'R', ('S', 1), pdf))

# sally_d = {'G': 0, 'H': 61, 'W': 110, 'S': 4, 'B': 22.78, 'T': 1, 'BR': 0.019, 'R': 0}
# print(cf.cf_closest_worlds(sally_d, ('S', 2), 'R'))

