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
            ('F', [], True),
            ('G', []),
			('A' , ['B', 'F']),
			('C', ['B', 'D']),
            ('D', ['A', 'G']),
            ('E', ['C'])
		]

varEquations = [
			    #'B = math.sin((t % 365) / 365 * 6.28) * 50 + 40 + normal(0, 5)',
                'B = exponential()',
                'F = logistic(5, 2)',
                'G = logistic(55, 10)',
			    'A = .5 * B + .3 * F + exponential()',
                'D = .25 * A + .35 * G + logistic(0,5)',
 			    'C = .1 * B + .2 * D + normal(0, 5)',
                'E = .5 * C + logistic(0,3)',
                't = t + 1'
		        ]
