import sys
import rv
import counterfactual
from synth import synthDataGen
from synth import getData

args = sys.argv

# If a test file is specified, only that file is run.
# Else we run the whole test suite specified in the tests var below.
if len(args) > 1:
    tests_synth = [args[1]]
else:
    tests_synth = ['cfBMI', 'cfHW']
    tests_synth = ['models/' + src + '.py' for src in tests_synth]
    tests_obs = []
    tests_obs = ['data/' + src for src in tests_obs]

bmi_sally = ('Sally', [('G', 0), ('H', 61), ('W', 110), ('S', 4), ('B', 22.78), ('T', 1), ('BR', 0.019), ('R', 0)])
hw_joe = ('Joe', [('X', 0.5), ('H', 1), ('Y', 1.5)])
queries = {'cfBMI': [(bmi_sally, 'R', ('S', 3)), (bmi_sally, 'R', ('S', 2)), (bmi_sally, 'R', ('S', 1))],
           'cfHW': [(hw_joe, 'Y', ('H', 2))]}

print('\n=============================')
print('Counterfactual Module Testing')
print('=============================\n')


def init_cf(test, obs):
    f = open(test, 'r')
    exec(f.read(), globals())

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

    # For data file, use the input file name with the .csv extension
    tokens = test.split('.')
    testFileRoot = str.join('.', tokens[:-1])
    datFileName = testFileRoot + '.csv'

    try:
        d = getData.DataReader(datFileName)
        data = d.read()
    except OSError as e:
        if obs:
            print('Error: Observational dataset not found.')
        else:
            print('Generating dataset...')
            synthDataGen.run(test, 100000)
            d = getData.DataReader(datFileName)
            data = d.read()

    return counterfactual.Counterfactual(gnodes, data, model, varEquations)


def run_tests(tests, desc, obs):
    for i, t in enumerate(tests):
        def run_test(test):
            cf = init_cf(test, obs)
            print('Now testing model:', test)
            print('Title:', testName)
            print('Description:', testDescript, '\n')
            for j, query in enumerate(queries[test[7:-3]]):
                print('\nQUERY #' + str(i+1) + '.' + str(j+1) + ': What would', query[0][0] + '\'s value of', query[1],
                      'be if', query[2][0], 'were', query[2][1], 'instead of', str(dict(query[0][1])[query[2][0]]) + '?')
                if 'HW' not in test:  # TODO remove when the hw model works
                    print('> Using the Data for Closest Worlds Computation:')
                    cf.closest_worlds(query[0][1], query[1], query[2], True)
                if model is not None and isLinear:
                    print('> Using the SEM for Deterministic Computation:')
                    cf.deterministic(query[0][1], query[1], query[2])
                # cf.probabilistic(query[0], query[1], query[2])
        if len(tests) > 1:
            print(desc, 'TEST #' + str(i+1))
            run_test(t)
            if i != len(tests) - 1:
                print('\n------------------------------------\n')
        else:
            run_test(t)


# Run the tests with synthetic data created from fully specified models
run_tests(tests_synth, 'SYNTHETIC', False)

# Run the tests with only observational data (model discovery TBD)
run_tests(tests_obs, 'REAL-WORLD', True)

# python counterfactualTest.py models/cfBMI.py
# python counterfactualTest.py models/cfHW.py

# from random import uniform
# import random
# random.seed(10)
# import math
#
# data = dict([('G', 0), ('H', 61), ('W', 110), ('S', 4), ('B', 22.78), ('T', 1), ('BR', 0.019), ('R', 0)])
#
# print(int((uniform(0,1) < (.9 -  data['BR'])/(1+math.log(data['S']))) if data['T'] else (uniform(0,1) <  (.5 - data['BR'])/(1+math.log(data['S'])))))
#
# # with cf:
# print(int((uniform(0,1) < (.9 -  data['BR'])/(1+math.log(2))) if data['T'] else (uniform(0,1) <  (.5 - data['BR'])/(1+math.log(2)))))