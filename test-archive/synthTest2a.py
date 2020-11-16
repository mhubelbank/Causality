# ----------------------------------------------------------------------------
# SEM Definitions
#
# Gaussian A->I set

varNames = ['A', 'B', 'C','D', 'E', 'F', 'G', 'H', 'I']
t=0
varEquations = [
				# 'A= normal(3, 5)',
				# 'B = normal(20,3)',
				# 'C = normal(-100,8)',
				# 'D = normal(10,5)',
				'A = gumbel(-45,7)',
				#'A = beta(.2,.7) * 5',
				#'B = math.sin(t)',
				'B = logistic(2,6)',
				'C = exponential() * 3',
				'D = exponential() * 4',
				'E = A - 1.2*B + 1.3 * C - .8 * D + logistic(-55,10)',
				'F = -1.1* E - logistic(100000, 10)',
				'G = .8*F + logistic(-10000,10)',
				'H = 1.2*F + logistic(0,10)',
				'I = -.9 * F + logistic(20,10)',
				't=t+.1',
				]
				
validation =  [('A', []), ('B', []), ('C', []), ('D', []), ('E', ['A','B','C','D']),
				('F', ['E']), ('G', ['F']), ('H',['F']), ('I', ['F']),
				]
# -----------------------------------------------------------------------------