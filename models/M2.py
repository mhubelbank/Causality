# ----------------------------------------------------------------------------
# Model Definition
# ----------------------------------------------------------------------------

# Initialize any variables here
t = 0

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
			    'B = noise()',
			    'A = coef() * B + noise()',
			    'C = coef() * A + coef() * B + noise()',
                't = t + 1'
		        ]
