# ----------------------------------------------------------------------------
# Model Definition
# ----------------------------------------------------------------------------

# Initialize any variables here
t = 0
from random import *
bRange = 2

# Describe the test
testDescript = 'Reference Model M2'

# Define the causal model.
# Each random variable has the following fields:
# - Name
# - List of Parents
# - isObserved (Optional, default True)
# - Data Type (Optional, default 'Numeric')
model =    [('B', []),
			('A' , ['B']),
			('C', ['B', 'A']),
			] 

# Structural Equation Model for data generation
varEquations = [
			    #'B = logistic(1, 1)',
			    #'A = .8 * B + normal(0, .5)',
			    #'C = 1.5 * A + 2 * B + normal(0, .2)',
                'B = logistic(3,5)',
			    'A = 1 * B  + logistic(0,3)',
			    'C = 1.5 * A + 2 * B + logistic(0,3)',
                't = t + 1'
		        ]
