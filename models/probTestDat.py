# ----------------------------------------------------------------------------
# Model Definition
# ----------------------------------------------------------------------------

# Data generator for testing the probability module (prob.py).
# Do not modify this file unless making corresponding changes to probTest.py
# Variables A - C are a deterministic emulation of two dice being rolled (ala Craps).
# A is die 1
# B is die 2
# C is the total roll:  die1 + die2
# N is a standardized normal distribution
# N2 is N + an offset normal distribution with mean = 1

# Initialize any variables here
t = 0
bRange = 2
# Describe the test
testDescript = 'Reference Data for probability testing.  Do not modify this file unless making corresponding changes to probTest.py'

# Define the causal model.
# Each random variable has the following fields:
# - Name
# - List of Parents
# - isObserved (Optional, default True)
# - Data Type (Optional, default 'Numeric')
model =    [('A', []),
			('B' , []),
            #('D', ['A']),
			('C', ['B', 'A']),
            ('N', []),
            ('N2', ['N'])
			]

# Structural Equation Model for data generation
varEquations = [
                #'B = choice(range(-bRange, bRange+1))',
                'A = choice(range(1, 7))',
			    #'B = 2 * B if uniform() < .2 else B',
                'B = choice(range(1, 7))',
                #'D = .5 * A + logistic(3,1)',
			    #'C = .5 * B + 1 * D + logistic(0,1)',
                'C = A + B',
                'N = normal(0,1)',
                'N2 = N + normal(1,1)',
                't = t + 1'
		        ]
