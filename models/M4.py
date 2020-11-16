# ----------------------------------------------------------------------------
# SEM Definitions
#

t = 0
testDescript = 'M4'

varNames = ['A', 'B', 'C', 'D']

varEquations = [
			    'B = math.sin((t % 365) / 365 * 6.28) * 50 + 40 + normal(0, 5)',
			    'A =  = .5 * B + normal(0,5)',
			    'C = .25 + A + .25 * B + + .2 * D + normal(0, 5)',
                'D = .25 * A + normal(0,5)'
                't = t + 1'
		        ]
				
model =    [('B', []),
			('A' , ['B']),
			('C', ['B', 'A', 'D']),
            ('D', ['A'])
			] 