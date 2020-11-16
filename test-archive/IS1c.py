# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Auto Regressive Integrated Moving Average (ARIMA(1,1,1)) Time-series model for Exogenous inputs with Hourglass Pattern Test (3 level)'

varNames = ['R11', 'R12',
			'R21',
			'R31', 'R32',
			]
R11_phi = [.7]
R11_theta = [.3]
R11_vals = [0] * 10
R11_errs = [0] * 10
R12_phi = [.7]
R12_theta = [.3]
R12_vals = [0] * 10
R12_errs = [0] * 10

varEquations = [
				'r11E = logistic(0,1)',
				'r12E = logistic(0,1)',
				'R11 = R11_vals[-1] + R11_phi[-1] * (R11_vals[-1] - R11_vals[-2]) + R11_theta[-1] * R11_errs[-1] + r11E',
				'R12 = R12_vals[-1] + R12_phi[-1] * (R12_vals[-1] - R12_vals[-2]) + R12_theta[-1] * R12_errs[-1] + r12E',
				'R21 = .5 * R11 + .5 * R12 + logistic(0,1)',
				'R31 = R21 + logistic(0,1)',
				'R32 = R21 + logistic(0,1)',
				'R11_vals.append(R11)',
				'temp = R11_vals.pop(0)',
				'R12_vals.append(R12)',
				'temp = R12_vals.pop(0)',
				'R11_errs.append(r11E)',
				'temp = R11_errs.pop(0)',
				'R12_errs.append(r12E)',
				'temp = R12_errs.pop(0)',
				]
				
validation =  [
				('R11', []), ('R12', []), ('R21', ['R11','R12']), ('R31', ['R21']), ('R32', ['R21']),
				]
# -----------------------------------------------------------------------------