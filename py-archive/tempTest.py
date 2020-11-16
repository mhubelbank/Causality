

from scipy.linalg import lstsq
import getData
import IndHSIC
import IndLn
import IndMI
import numpy as np
import synthDataGen
from debug import *
import math
import time


#IND_TYPE = 'HSIC'
IND_TYPE = 'PL' # Pairwise Lingam
#IND_TYPE = 'Ln'
#IND_TYPE = 'MI'
firstTime = True

RUNS = 100
DATA_POINTS = 1000
MAX_DIFFICULTY = 0
FILE = 'synthTest4'

if IND_TYPE == 'PL':
	import CaLingam
	cal = CaLingam.IC('..\Tests\\' + FILE + '.csv', 1000)

SNR = 10**100
SNR_DB = 0
CUM_TIME = 0

def run(s1, s2):
	global firstTime
	global SNR, SNR_DB, CUM_TIME
	reset = True
	if firstTime:
		#print('std ratio = ', s1D.std() / s2D.std())
		firstTime = False
	else:
		reset = False
	synthDataGen.run('..\Tests\\' + FILE + '.py', samples=DATA_POINTS, maxDifficulty=MAX_DIFFICULTY, reset=reset)
	dr = getData.DataReader(input='..\Tests\\' + FILE + '.csv', limit=DATA_POINTS)
	s1D = np.array(dr.getSeries(s1))
	s2D = np.array(dr.getSeries(s2))
	#print('s1 = ', s1D.mean(), s1D.std())
	#print('s2 = ', s2D.mean(), s2D.std())
	if firstTime:
		#print('std ratio = ', s1D.std() / s2D.std())
		firstTime = False
	coefs = np.polyfit(s1D, s2D.T, 1)
	x1coef, x0coef = coefs
	Dprint('coefs = ', coefs)
	cs2D = s1D * x1coef + x0coef

	# coef = lstsq(np.array([list(s1D)]).T, s2D)[0][0]
	# Dprint('coef = ', coef)
	# cs2D = s1D * coef

	res = s2D - cs2D
	signal = cs2D.var(ddof=1)
	noise = res.var(ddof=1)
	snr = signal / noise
	if snr < SNR:
		SNR = snr
		SNR_DB = 10 * math.log(snr,10)

	start = time.time()
	if IND_TYPE == 'PL':
		dep = scoreDependence(s1D, s2D)
	else:
		dep = scoreDependence(s1D, res)
	end = time.time()
	duration = end - start
	CUM_TIME += duration
	Dprint('Residual Dependence for ', s1, '-->', s2, ' = ', dep)
	return dep

def scoreDependence(s1, s2):
	if IND_TYPE == 'HSIC':
		return IndHSIC.scoreDependence(s1, s2)
	elif IND_TYPE == 'Ln':
		return IndLn.scoreDependence(s1, s2)
	elif IND_TYPE == 'PL':
		dir = cal.direction2(s1, s2)
		return 1 - dir
	elif IND_TYPE == 'MI':
		return IndMI.scoreDependence(s1, s2)
		#print('dir = ', dir)
		
import sys
s1 = sys.argv[1]
s2 = sys.argv[2]

print('FILE = ', FILE)
print('RUNS = ', RUNS)
print('IND_TYPE = ', IND_TYPE)
print('DATA_POINTS = ', DATA_POINTS)

scores = {s1:0, s2:0}
#IND_TYPE = method
for i in range(RUNS):
	dep1 = run(s1, s2)
	dep2 = run(s2, s1)
	if dep1 < dep2:
		Dprint('predecessor = ', s1)
		scores[s1] = scores[s1] + 1
	else:
		Dprint('predecessor = ', s2)
		scores[s2] = scores[s2] + 1
sc1 = scores[s1]
sc2 = scores[s2]

Dprint('FINAL SCORES = ', sc1, ' / ', sc2)
print('DIFFICULTY = ', SNR_DB)
normDiff = SNR_DB/(DATA_POINTS**.5) * 100
print('NORM DIFFICULTY = ', normDiff)
if sc1 > sc2:
	print('predecessor = ', s1, '(', sc1/(sc1 + sc2)*100, '%)')
else:
	print('predecessor = ', s2, '(', sc2/(sc1 + sc2)*100, '%)')
power = max([((sc1 / float(RUNS) * 100) -50) / 50.0, 0.0])
print('POWER = ', power)
if power >= 0.1 and power <= .9 and normDiff > 0.0:
	scaledPower = (((math.e**power)-1) / (math.e - 1))**.5
	print('EFFECTIVENESS = ', scaledPower * SNR_DB)
	sensitivity = scaledPower * normDiff
	print('SENSITIVITY = ', sensitivity)
	print('RUN TIME = ', CUM_TIME)
	print('DESIRABILITY = ', sensitivity / CUM_TIME)
	#methScores[method] = ((sc1 / float(RUNS) * 100.0) - 50.0)/50.0
	#scHSIC = methScores['HSIC']
	#print('methScores = ', methScores)
	#scLn = methScores['Ln']
	# if scHSIC > scLn:
		# winner = scHSIC
		# loser = scLn
		# winName = 'HSIC'
	# else:
		# winner = scLn
		# loser = scHSIC
		# winName = 'Ln'
	# print('Best algorithm = ', winName, '(', winner / loser, ')')
else:
	print('Other metrics cannot be computed with power less than .1 or power > .9 or negative Difficulty')
	
	
