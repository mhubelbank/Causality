# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Moving Average Time-series model for Exogenous inputs with Hourglass Pattern Test (3 level)'

varNames = ['R11', 'R12',
			'R21',
			'R31', 'R32',
			]

lastR11E = 0.0
theta11 = .7
lastR12E = 0.0
theta12 = .3

varEquations = [
				'r11E = logistic(0,1)',
				'r12E = logistic(0,1)',
				'R11 = theta11 * lastR11E + r11E',
				'R12 = theta12 * lastR12E + r12E',
				'R21 = .5 * R11 + .5 * R12 + logistic(0,1)',
				'R31 = R21 + logistic(0,1)',
				'R32 = R21 + logistic(0,1)',
				'lastR11E = r11E',
				'lastR12E = r12E',
				]
				
validation =  [
				('R11', []), ('R12', []), ('R21', ['R11','R12']), ('R31', ['R21']), ('R32', ['R21']),
				]
# -----------------------------------------------------------------------------