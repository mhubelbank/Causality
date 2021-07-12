# ----------------------------------------------------------------------------
# Model Definition
# ----------------------------------------------------------------------------

# Initialize any variables here
# bRange = 2
# Describe the test
testDescript = 'Counterfactual Test Model Det'

# Define the causal model.
# Each random variable has the following fields:
# - Name
# - List of Parents
# - isObserved (Optional, default True)
# - Data Type (Optional, default 'Numeric')
model = [('X', []),
         ('H', ['X']),
         ('Y', ['X', 'H']),
         ]

# Structural Equation Model for data generation
varEquations = [
    'X = normal(0, 5)',
    'H = 0.5 * X',
    'Y = 0.7 * X + 0.4 * H'
]

# python synth/synthDataGen.py models/cfTestDet.py 100000