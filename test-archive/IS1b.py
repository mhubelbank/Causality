# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Auto Regressive Moving Average (ARMA) Time-series model for Exogenous inputs with Hourglass Pattern Test (3 level)'

varNames = ['R11', 'R12',
			'R21',
			'R31', 'R32',
			]

lastR11 = 0.0
lastR12 = 0.0
phi11 = .7
phi12 = .3
theta11 = .3
theta12 = .7
lastR11E = 0.0
lastR12E = 0.0


varEquations = [
				'r11E = logistic(0,1)',
				'r12E = logistic(0,1)',
				'R11 = phi11 * lastR11 + theta11 * lastR11E + r11E',
				'R12 = phi12 * lastR12 + theta12 * lastR12E + r12E',
				'R21 = .5 * R11 + .5 * R12 + logistic(0,1)',
				'R31 = R21 + logistic(0,1)',
				'R32 = R21 + logistic(0,1)',
				'lastR11E = r11E',
				'lastR12E = r12E',
				'lastR11 = R11',
				'lastR12 = R12',
				]
				
validation =  [
				('R11', []), ('R12', []), ('R21', ['R11','R12']), ('R31', ['R21']), ('R32', ['R21']),
				]
# -----------------------------------------------------------------------------