# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Disjoint Structures Test (V-Structure, IV-Structure, Chain, Unconnected)'

varNames = ['A1', 'A2', 'A3', 
			'B1', 'B2', 'B3',
			'C1', 'C2', 'C3',
			'D1', 'D2', 'D3',
			]

varEquations = [
				# V-Structure
				'A1 = logistic(0,1)',
				'A2 = logistic(0,1)',
				'A3 = .5 * A1 + .5 * A2 + logistic(0,1)',
				# IV-Structure
				'B1 = logistic(0,1)',
				'B2 = B1 + logistic(0,1)',
				'B3 = B1 + logistic(0,1)',
				# Chain
				'C1 = logistic(0,1)',
				'C2 = C1 + logistic(0,1)',
				'C3 = C2 + logistic(0,1)',
				# Unconnected
				'D1 = logistic(0,1)',
				'D2 = logistic(0,1)',
				'D3 = logistic(0,1)',
				]
				
validation =  [
				('A1', []), ('A2', []), ('A3', ['A1','A2']), 
				('B1', []), ('B2', ['B1']), ('B3', ['B1']),
				('C1', []), ('C2',  ['C1']), ('C3', ['C2']),
				('D1', []), ('D2', []), ('D3', []),
				]
# -----------------------------------------------------------------------------