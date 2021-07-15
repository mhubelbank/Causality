# ----------------------------------------------------------------------------
# Model Definition
# ----------------------------------------------------------------------------

# Initialize any variables here
# bRange = 2
# Describe the test
testName = 'BMI, Treatment Stage on Chances of Disease Recovery'
isLinear = False
testDescript = ''

# Define the causal model.
# Each random variable has the following fields:
# - Name
# - List of Parents
# - isObserved (Optional, default True)
# - Data Type (Optional, default 'Numeric')
model = [('G', []),
         ('H', ['G']),
         ('W', ['G']),
         ('S', []),
         ('B', ['W', 'H']),
         ('T', []),
         ('BR', ['B']),
         ('R', ['BR', 'S']),
         ]

# Structural Equation Model for data generation
varEquations = [
    'G = choice([0,1])',
    'H = 65 + normal(0, 3) + G * 6',
    'W = 100 + abs(normal(30, 30))  + G * 30',
    'S = choice([1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4])',
    'B = (W / (H**2) * 703 + normal(0, 4))',
    'T = choice([0,1])',
    'BR = abs((22.5 - B)/22.5) * 1.5',
    'R = int((uniform(0,1) < (.9 -  BR)/(1+math.log(S))) if T else (uniform(0,1) <  (.5 - BR)/(1+math.log(S))))'
]
# If perfect weight, treated, stage 1, 90% chance of recovery
# IF perfect weight, untreated, stage 1, 50% chance of recovery
# Stage and BMI affect recovery chances; stage: higher = better, BMI: too high or too low decreases chances.

# python synth/synthDataGen.py models/cfTestBMI.py 100000