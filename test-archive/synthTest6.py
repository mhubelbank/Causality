# ----------------------------------------------------------------------------
# Simplified version of Test1
# SEM Definitions

varNames = ['A', 'B', 'C','D', 'E', 'F', 'G', 'H', 'I']

varEquations = ['A=  data()',
				'B =  data()',
				'C = data()',
				'D = data()',
				'E = coef() * A + coef() * B + coef() * C + coef() * D + noise()',
				'F = coef() * E + noise()',
				'G = coef() * F + noise()',
				'H = coef() * F + noise()',
				'I = coef() * F + noise()',
				]
				
validation =  [('A', []), ('B', []), ('C', []), ('D', []), ('E', ['A','B','C','D']),
				('F', ['E']), ('G', ['F']), ('H',['F']), ('I', ['F']), 
				]

# -----------------------------------------------------------------------------