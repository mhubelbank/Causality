# ----------------------------------------------------------------------------
# Model Definition
# ----------------------------------------------------------------------------

# Initialize any variables here
# bRange = 2
# Describe the test
testName = 'Homework on Exam Score'
isLinear = True
testDescript = ''

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
varEquations = [  # Each value given as the number of standard deviations above the mean the student falls
    'X = normal(0, 5)',  # Time spent in an after-school remedial program
    'H = 0.5 * X',  # Amount of homework completed
    'Y = 0.7 * X + 0.4 * H'  # Score on exam
]

# python synth/synthDataGen.py models/cfHW.py 100000