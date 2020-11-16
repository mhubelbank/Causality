# ----------------------------------------------------------------------------
# SEM Definitions
#

t = 0
testDescript = 'Reference Model M5'

varNames = ['A', 'B', 'C', 'D', 'E']

varEquations = [
			    'B = math.sin((t % 365) / 365 * 6.28) * 50 + 40 + normal(0, 5)',
			    'A =  = .5 * B + normal(0,5)',
			    'C = .25 + A + .25 * B + + .2 * D + normal(0, 5)',
                'D = .25 * A + normal(0,5)',
                'E = .5 * C + normal(0,3)',
                't = t + 1'
		        ]
				
model = [   ('B', []),
			('A' , ['B']),
			('C', ['B', 'D']),
            ('D', ['A']),
            ('E', ['C'])
		]