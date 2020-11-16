# ----------------------------------------------------------------------------
# SEM Definitions

varNames = ['ZA', 'B', 'C','D', 'AE', 'F', 'G', 'H', 'DI']

varEquations = ['ZA= logistic(0,1)',
				'B = logistic(0,1)',
				'C = logistic(0,1)',
				'D = logistic(0,1)',
				'AE = ZA + B + C + D + logistic(0,1)',
				'F = AE + logistic(0,1)',
				'G = F + logistic(0,1)',
				'H = F + logistic(0,1)',
				'DI = F + logistic(0,1)',
				]
				
validation =  [('ZA', []), ('B', []), ('C', []), ('D', []), ('AE', ['ZA','B','C','D']),
				('F', ['AE']), ('G', ['F']), ('H',['F']), ('DI', ['F']), 
				]
# -----------------------------------------------------------------------------