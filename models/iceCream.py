# ----------------------------------------------------------------------------
# SEM Definitions
#

t = 0
testDescript = 'Ice Cream Test'

varNames = ['temperature', 'iceCream', 'crime']

varEquations = [
			    'temperature = math.sin((t % 365) / 365 * 6.28) * 50 + 40 + normal(0, 5)',
			    'iceCream = 100000 + 10 * temperature + normal(0, 1000)',
			    'crime = 10 + temperature + normal(0, 10)',
                't = t + 1'
		        ]
				
model =    [('temperature', []),
			('iceCream' , ['temperature']),
			('crime', ['temperature'])
			] 