# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Auto-regression (AR) Time-Series Model -- Hourglass Pattern Test (3 level)'

varNames = ['R11', 'R12',
			'R21',
			'R31', 'R32',
			]

lastR11 = 0.0
phi11 = .7
lastR12 = 0.0
phi12 = .3

varEquations = [
				'R11 = phi11 * lastR11 + logistic(0,1)',
				'R12 = phi12 * lastR12 + logistic(0,1)',
				'R21 = .5 * R11 + .5 * R12 + logistic(0,1)',
				'R31 = R21 + logistic(0,1)',
				'R32 = R21 + logistic(0,1)',
				'lastR11 = R11',
				'lastR12 = R12',
				]
				
validation =  [
				('R11', []), ('R12', []), ('R21', ['R11','R12']), ('R31', ['R21']), ('R32', ['R21']),
				]
# -----------------------------------------------------------------------------