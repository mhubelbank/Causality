# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Auto-Regressive (AR) Hourglass Pattern Test (3 level)'

varNames = ['R11', 'R12',
			'R21',
			'R31', 'R32',
			]

phi11 = .7
prevR11 = 0.0
phi12 = .3
prevR12 = 0.0

varEquations = [
				'R11 = phi11 * prevR11 + normal(0,1)',
				'R12 = phi12 * prevR12 + normal(0,1)',
				'R21 = .5 * R11 + .5 * R12 + normal(0,1)',
				'R31 = R21 + normal(0,1)',
				'R32 = R21 + normal(0,1)',
				'prevR11 = R11',
				'prevR12 = R12',
				]
				
validation =  [
				('R11', []), ('R12', []), ('R21', ['R11','R12']), ('R31', ['R21']), ('R32', ['R21']),
				]
# -----------------------------------------------------------------------------