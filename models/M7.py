# ----------------------------------------------------------------------------
# SEM Definitions
#
# ----------------------------------------------------------------------------
# Model Definition
# ----------------------------------------------------------------------------

# Initialize any variables here
t = 0

# Describe the test
testDescript = 'Reference Model M7'

# Define the causal model.
# Each random variable has the following fields:
# - Name
# - List of Parents
# - isObserved (Optional, default True)
# - Data Type (Optional, default 'Numeric')  
model = [   ('B', []),
            ('F', [], False),
            ('G', []),
			('A' , ['B', 'F']),
			('C', ['B', 'D']),
            ('D', ['A', 'G']),
            ('E', ['C'])
		]

varEquations = [
			    'B = math.sin((t % 365) / 365 * 6.28) * 50 + 40 + normal(0, 5)',
                'F = normal(10, 7)',
                'G = logistic(55, 10)',
			    'A = .5 * B + .3 * F + normal(0,5)',
                'D = .25 * A + .35 * G + normal(0,5)',
 			    'C = .25 * B + .2 * D + normal(0, 5)',
                'E = .5 * C + normal(0,3)',
                't = t + 1'
		        ]
