# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Auto Regressive Integrated Moving Average (ARIMA(1,1,1)) Time-series model for Exogenous inputs with Hourglass Pattern Test (3 level)'

varNames = ['R11',
			'R21',
			]
R11_phi = [.7]
R11_theta = [.3]
R11_vals = [0] * 10
R11_errs = [0] * 10


varEquations = [
				'r11E = logistic(0,1)',
				'R11 = R11_vals[-1] + R11_phi[-1] * (R11_vals[-1] - R11_vals[-2]) + R11_theta[-1] * R11_errs[-1] + r11E',
				#'R11 = R11_phi[-1] * R11_vals[-1] + R11_theta[-1] * R11_errs[-1] + r11E',
				'R21 = R11 + logistic(0,1)',
				'R11_vals.append(R11)',
				'temp = R11_vals.pop(0)',
				'R11_errs.append(r11E)',
				'temp = R11_errs.pop(0)',
				]
				
validation =  [
				('R11', []), ('R21', ['R11']), 
				]
# -----------------------------------------------------------------------------