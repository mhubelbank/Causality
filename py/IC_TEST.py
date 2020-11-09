# IC_TEST
#
# Master Test for the Inferred Causality (IC) Prototype Subsystem

# PARAMS are: (<datafile>, <runs>, <validationFunction>)
from debug import *
import CaLingam
import CaPC
import CaCCM
import CaICAny
import synthDataGen
import sys
import time
import analyzeSEM
import fnmatch

VALIDATION_ERRORS = {}
CAUSAL_ORDERS = {}

TESTS = [
		# TestDef := (<testName>, <testType := 'synth' | 'live'>, <algorithm (see below)> <testDescription>, <regenCount or 0>, <maxDifficulty or 0>, [<Params>])
		# regenCount :=  Number of times to run before resetting the data generator to create a new parameterization.  Zero to disable.
		# maxDifficulty := The maximum difficulty for the data generator to generate or zero to disable.  Note: This can cause a reduction in
		#    performance of the tests, or even stalling if attempting to create tests with too small a maxDifficulty.
		# algorithm := 'PC' (PC independence algorithm ala Pearl), 'LI' (LiNGAM algo based on Darmois-Skitovich Theorem with default method), 
		#	'MI' (LiNGAM with Mutual Information independence method), 'PL' (Pairwise Lingam),  
		#   'CM' (CauseMap non-linear CCM based on Taken's Theorem), 'ANY' (Choose the best method automatically)

		# New Rudimentals Tests
		('R1', 'synth', 'MI', 'V-Structure', 0, 0, [
			('../Tests/R1.py', 10, 10000, 'SynthStructVal'),
			]),
		('R1a', 'synth', 'ANY', 'V-Structure 3->1', 0, 0, [
			('../Tests/R1a.py', 10, 10000, 'SynthStructVal'),
			]),
		('R1b', 'synth', 'LI', 'V-Structure 5->1', 0, 0, [
			('../Tests/R1b.py', 10, 10000, 'SynthStructVal'),
			]),
		('R2', 'synth', 'LI', 'Inverted V-Structure', 0, 0, [
			('../Tests/R2.py', 10, 10000, 'SynthStructVal'),
			]),
		('R2a', 'synth', 'LI', 'Inverted V-Structure 1->3', 0, 0, [
			('../Tests/R2a.py', 10, 10000, 'SynthStructVal'),
			]),
		('R2b', 'synth', 'LI', 'Inverted V-Structure 1->5', 0, 0, [
			('../Tests/R2b.py', 10, 10000, 'SynthStructVal'),
			]),
		('R3', 'synth', 'LI', '3-Chain', 0, 0, [
			('../Tests/R3.py', 10, 10000, 'SynthStructVal'),
			]),
		('R3a', 'synth', 'LI', '4-Chain', 0, 0, [
			('../Tests/R3a.py', 10, 10000, 'SynthStructVal'),
			]),
		('R3b', 'synth', 'LI', '7-Chain', 0, 0, [
			('../Tests/R3b.py', 10, 10000, 'SynthStructVal'),
			]),
		('R4', 'synth', 'LI', 'W-Pattern Test', 0, 0, [
			('../Tests/R4.py', 10, 10000, 'SynthStructVal'),
			]),
		('R5', 'synth', 'LI', 'Inverted W-Pattern Test', 0, 0, [
			('../Tests/R5.py', 10, 10000, 'SynthStructVal'),
			]),
		('R6', 'synth', 'LI', 'Triangle Pattern Test (3 level)', 0, 0, [
			('../Tests/R6.py', 10, 10000, 'SynthStructVal'),
			]),
		('R6a', 'synth', 'LI', 'Triangle Pattern Test (4 level)', 0, 0, [
			('../Tests/R6a.py', 10, 10000, 'SynthStructVal'),
			]),
		('R7', 'synth', 'LI', 'Inverted Triangle Pattern Test (3 level)', 0, 0, [
			('../Tests/R7.py', 10, 10000, 'SynthStructVal'),
			]),
		('R7a', 'synth', 'LI', 'Inverted Triangle Pattern Test (4 level)', 0, 0, [
			('../Tests/R7a.py', 10, 10000, 'SynthStructVal'),
			]),
		('R8', 'synth', 'LI', 'Diamond Pattern Test (3 level)', 0, 0, [
			('../Tests/R8.py', 10, 10000, 'SynthStructVal'),
			]),
		('R8a', 'synth', 'LI', 'Diamond Pattern Test (5 level)', 0, 0, [
			('../Tests/R8a.py', 10, 10000, 'SynthStructVal'),
			]),
		('R8b', 'synth', 'LI', 'Diamond Pattern Test (7 level)', 0, 0, [
			('../Tests/R8b.py', 1, 10000, 'SynthStructVal'),
			]),
		('R9', 'synth', 'LI', 'Hourglass Pattern Test (3 level)', 0, 0, [
			('../Tests/R9.py', 10, 10000, 'SynthStructVal'),
			]),
		('R9a', 'synth', 'PC', 'Hourglass Pattern Test (5 level)', 0, 0, [
			('../Tests/R9a.py', 10, 10000, 'SynthStructVal'),
			]),
		('R10', 'synth', 'LI', 'Cascade Pattern Test (3 level)', 0, 0, [
			('../Tests/R10.py', 10, 10000, 'SynthStructVal'),
			]),
		('R10a', 'synth', 'LI', 'Cascade Pattern Test (4 level)', 0, 0, [
			('../Tests/R10a.py', 10, 10000, 'SynthStructVal'),
			]),
		('R11', 'synth', 'LI', 'Disjoint Structure Test (V, IV, Chain, Unconnected)', 0, 0, [
			('../Tests/R11.py', 10, 10000, 'SynthStructVal'),
			]),		
		# Rudimentals Gaussian Variants
		('RG1', 'synth', 'PC', 'V-Structure (Gaussian)', 0, 0, [
			('../Tests/RG1.py', 10, 10000, 'SynthStructVal'),
			]),
		('RG2', 'synth', 'PC', 'Inverted V-Structure (Gaussian)', 0, 0, [
			('../Tests/RG2.py', 10, 100000, 'SynthStructVal'),
			]),		
		('RG3', 'synth', 'PC', '3-Chain (Gaussian)', 0, 0, [
			('../Tests/RG3.py', 10, 100000, 'SynthStructVal'),
			]),
		('RG4', 'synth', 'PC', 'Diamond Pattern Test (Gaussian)', 0, 0, [
			('../Tests/RG4.py', 10, 10000, 'SynthStructVal'),
			]),
		('RG9', 'synth', 'PC', 'Hourglass Pattern Test (3 level)', 0, 0, [
			('../Tests/RG9.py', 10, 10000, 'SynthStructVal'),
			]),
		('RG9a', 'synth', 'PC', 'Hourglass Pattern Test (5 level)', 0, 0, [
			('../Tests/RG9a.py', 10, 10000, 'SynthStructVal'),
			]),

		# Rudimentals Non-Linear Variants
		('RN9', 'synth', 'CM', 'Hourglass Pattern Test (3 level)', 0, 0, [
			('../Tests/RN9.py', 10, 100000, 'SynthStructVal'),
			]),

		# Input Sensitivity Tests
		('IS1', 'synth', 'LI', 'Auto-Regressive (AR(1)) exogenous inputs Hourglass Structure', 0, 0, [
			('../Tests/IS1.py', 10, 10000, 'SynthStructVal'),
			]),
		('IS1a', 'synth', 'LI', 'Moving Average (MA(1)) exogenous inputs Hourglass Structure', 0, 0, [
			('../Tests/IS1a.py', 10, 10000, 'SynthStructVal'),
			]),
		('IS1b', 'synth', 'LI', 'Auto-Regressive Moving Average (ARMA(1,1)) exogenous inputs Hourglass Structure', 0, 0, [
			('../Tests/IS1b.py', 10, 10000, 'SynthStructVal'),
			]),
		('IS1c', 'synth', 'LI', 'Auto-Regressive Integrated Moving Average (ARIMA(1,1,1)) exogenous inputs Hourglass Structure', 0, 0, [
			('../Tests/IS1c.py', 10, 10000, 'SynthStructVal'),
			]),
		('IS1d', 'synth', 'LI', 'Auto-Regressive Integrated Moving Average (ARIMA(1,1,1)) exogenous inputs Hourglass Structure', 0, 0, [
			('../Tests/IS1d.py', 100, 10000, 'SynthStructVal'),
			]),
		# Input Sensitivity Tests (Gaussian Variants)
		('ISG1', 'synth', 'PC', 'Auto-Regressive (AR) exogenous inputs Hourglass Structure', 0, 0, [
			('../Tests/ISG1.py', 10, 100000, 'SynthStructVal'),
			]),
		('ISG1a', 'synth', 'PC', 'Moving Average (MA) exogenous inputs Hourglass Structure', 0, 0, [
			('../Tests/ISG1a.py', 10, 100000, 'SynthStructVal'),
			]),
		('ISG1b', 'synth', 'PC', 'Auto-Regressive Moving Average (ARMA) exogenous inputs Hourglass Structure', 0, 0, [
			('../Tests/ISG1b.py', 10, 100000, 'SynthStructVal'),
			]),
		('IS2', 'synth', 'LI', 'Random walk exogenous single input to IV Structure', 0, 0, [
			('../Tests/IS2.py', 10, 10000, 'SynthStructVal'),
			]),
		('IS2G', 'synth', 'PC', 'Random walk exogenous single input to IV Structure (Gaussian)', 0, 0, [
			('../Tests/IS2G.py', 10, 100000, 'SynthStructVal'),
			]),
		('IS3', 'synth', 'LI', '3-Chain with sine wave for exogenous variable', 0, 0, [
			('../Tests/IS3.py', 10, 100000, 'SynthStructVal'),
			]),
		# Transformation Sensitivity Tests
		
		# Legacy Tests
		('Test1', 'synth', 'LI', 'Basic Non-Gaussian A->I Pattern', 0, 0, [
			# Params := (<DataGeneratorFile(synth) or DataFile(live)>, <how many times to run>, <how many data points per variable>, <validation function>, <validation data override(optional-->uses validation from input file by default)

			('../Tests/synthTest1.py', 10, 10000, 'SynthStructVal'),
			('../Tests/synthTest1.py', 10, 10000, 'SynthOrderVal'),
			('../Tests/synthTest1.py', 100, 5000, 'SynthOrderVal'),
			('../Tests/synthTest1.py', 100, 1000, 'SynthOrderVal'),
			('../Tests/synthTest1.py', 100, 500, 'SynthOrderVal'),
			('../Tests/synthTest1.py', 100, 100, 'SynthOrderVal'),
			]),
		('Test2', 'synth', 'ANY', 'Gaussian A->I Pattern', 0, 0, [
			('../Tests/synthTest2.py', 10, 10000, 'SynthStructVal'),
			('../Tests/synthTest2.py', 10, 10000, 'SynthOrderVal'),
			('../Tests/synthTest2.py', 100, 100000, 'SynthStructVal'),
			('../Tests/synthTest2.py', 10, 50000, 'SynthOrderVal'),
			('../Tests/synthTest2.py', 10, 50000, 'SynthStructVal'),
			('../Tests/synthTest2.py', 100, 10000, 'SynthStructVal'),
			('../Tests/synthTest2.py', 100, 5000, 'SynthStructVal'),
			('../Tests/synthTest2.py', 100, 1000, 'SynthStructVal'),
			]),
		('Test2a', 'synth', 'ANY', 'A -> Pattern using similar parameters to Test 2, but non-gaussian, using LiNGAM', 0,0,[
			('../Tests/synthTest2a.py', 10, 10000, 'SynthStructVal'),
			]),
		('Test2b', 'synth', 'LI', 'Gaussian A->I Pattern using LiNGAM (Should fail)', 0, 0, [
			('../Tests/synthTest2.py', 10, 10000, 'SynthStructVal'),
			]),
		('Test2c', 'synth', 'PC', 'Non-Gaussian A->I Pattern using PC', 0, 0, [
			('../Tests/synthTest2a.py', 10, 100000, 'SynthStructVal'),
			('../Tests/synthTest2a.py', 10, 50000, 'SynthStructVal'),
			('../Tests/synthTest2a.py', 10, 25000, 'SynthStructVal'),
			('../Tests/synthTest2a.py', 10, 10000, 'SynthStructVal'),
			('../Tests/synthTest2a.py', 10, 5000, 'SynthStructVal'),
			('../Tests/synthTest2a.py', 10, 1000, 'SynthStructVal'),
			]),
		

		('Test3', 'synth', 'CM', 'Basic A->I with non-linear terms', 0, 0, [
			('../Tests/synthTest3.py', 10, 100000, 'SynthOrderVal'),
			]),
		('Test4', 'synth', 'LI', 'Basic A->I with simulated time-series i.e., non-iid', 0, 0, [
			('../Tests/synthTest4.py', 100, 10000, 'SynthOrderVal'),
			('../Tests/synthTest4.py', 100, 500, 'SynthStructVal'),
			('../Tests/synthTest4.py', 100, 100, 'SynthStructVal'),
			('../Tests/synthTest4.py', 100, 50, 'SynthStructVal'),
			('../Tests/synthTest4.py', 100, 10, 'SynthStructVal'),
			]),
		('Test5', 'synth', 'LI', 'Experimental', 0, 0, [
			('../Tests/synthTest5.py', 20, 1000, 'SynthStructVal'),
			]),
			
		('Test6', 'synth', 'LI', 'A->I with random parameterization', 10, 2.5, [
			('../Tests/synthTest6.py', 100, 1000, 'SynthOrderVal'),
			('../Tests/synthTest6.py', 100, 8000, 'SynthOrderVal'),
			('../Tests/synthTest6.py', 100, 5000, 'SynthOrderVal'),
			('../Tests/synthTest6.py', 100, 3000, 'SynthOrderVal'),
			('../Tests/synthTest6.py', 100, 1000, 'SynthOrderVal'),		
			]),
		('Test7', 'synth', 'LI', 'Two variable test with random parameterization', 10, 25.0, [
			('../Tests/synthTest7.py', 100, 100, 'SynthOrderVal'),
			('../Tests/synthTest7.py', 100, 200, 'SynthOrderVal'),
			('../Tests/synthTest7.py', 100, 500, 'SynthOrderVal'),
			('../Tests/synthTest7.py', 100, 1000, 'SynthOrderVal'),
			('../Tests/synthTest7.py', 100, 2000, 'SynthOrderVal'),
			('../Tests/synthTest7.py', 100, 5000, 'SynthOrderVal'),
			]),
		('Test8', 'synth', 'LI', 'Dougs A->I test with all coefs and stddevs = 1', 0, 0, [
			('../Tests/synthTest8.py', 10, 10000, 'SynthOrderVal'),
			('../Tests/synthTest8.py', 10, 8000, 'SynthOrderVal'),
			('../Tests/synthTest8.py', 10, 5000, 'SynthOrderVal'),
			('../Tests/synthTest8.py', 10, 3000, 'SynthOrderVal'),
			('../Tests/synthTest8.py', 10, 1000, 'SynthOrderVal'),
			]),
		('Test9', 'synth', 'LI', 'Experimental', 0, 0, [
			('../Tests/synthTest9.py', 100, 20000, 'SynthOrderVal'),
			('../Tests/synthTest9.py', 200, 10000, 'SynthOrderVal'),
			('../Tests/synthTest9.py', 400, 5000, 'SynthOrderVal'),
			('../Tests/synthTest9.py', 600, 2000, 'SynthOrderVal'),
			('../Tests/synthTest9.py', 800, 1000, 'SynthOrderVal'),
			('../Tests/synthTest9.py', 1000, 100, 'SynthOrderVal'),
			]),
		('Test10', 'synth', 'LI', 'Non-alphabetic order test', 0, 0, [
			('../Tests/synthTest10.py', 10, 10000, 'SynthOrderVal'),
			]),
		('Test12', 'synth', 'LI', 'Non-alphabetic order test', 0, 0, [
			('../Tests/synthTest12.py', 10, 10000, 'SynthStructVal'),
			]),
		('Sim1', 'live', 'LI', 'Data from Dougs Network Simulator', 0, 0, [
			('../Tests/simulation1.csv', 100, 1000, ''),
			]),
		('Ex1', 'live', 'LI', 'Live data from Mikes network with test points A-J', 0, 0, [
			('../data/experiment1.csv', 20, 2000, ''),
			]),
		]
		
def isTestSelected(selectionList, testName):
	if selectionList is not None:
		for term in selectionList:
			if fnmatch.fnmatch(testName, term):
				return True
	return False
def runTests(selectionList):
	global RESET_COUNT, MAX_DIFFICULTY
	for test in TESTS:
		testName, testType, testAlgo, testDescr, resetCount, maxDifficulty, testParams = test[:7]
		RESET_COUNT = resetCount
		MAX_DIFFICULTY = maxDifficulty
		if selectionList is None or len(selectionList) == 0 or isTestSelected(selectionList, testName):
			print('Starting Test: ', testName, '(', testDescr, ')', '.  Parms = ', test[:6])
			runTest(testName, testType, testAlgo, testDescr, testParams)
			print('Test: ', testName, ' Complete/n')
	return
	
def runTest(testName, testType, testAlgo, testDescr, testParams):
		
	for paramSet in testParams:
		# Run through one test
		#print()
		print('Running test with parms = ', paramSet)
		datafile, runs, validation = paramSet[:3]
		CaLingamTest(testType, testAlgo, paramSet)
	return

	
RESET_COUNT = 0
MAX_DIFFICULTY = 0

def CaLingamTest(testType, testAlgo, paramSet):
	global VALIDATION_ERRORS, CAUSAL_ORDERS
	VALIDATION_ERRORS = {}
	CAUSAL_ORDERS = {}
	DAGS = {}
	maxDifficulty = MAX_DIFFICULTY
	fails = 0.0
	datafile, runs, dataPoints, validation = paramSet[:4]
	valData = None
	totalDuration = 0
	datafileRootName = datafile.split('.')[-2].split('//')[-1]
	#print('dfrn = ', datafileRootName)
	sa = analyzeSEM.SemAnalyzer(datafile)
	difficulty = 0
	diffNorm = 0
	for i in range(runs):
		if i == 0:
			reset = True
		if RESET_COUNT > 0 and i > 0 and i / RESET_COUNT == int(i / RESET_COUNT):
			reset = True
			print()
			print('Previous SEM:')
			print(synthDataGen.getSEM())
			#print('Resetting')
		else:
			reset = False
		if i == 0:
			reset = True
		prune = True
		if validation == 'SynthOrderVal':
			# Suppress pruning in order to increase perf when only testing order.
			prune = False
		if testType == 'synth':
			outFile = synthDataGen.run(datafile, samples=dataPoints, reset=reset, maxDifficulty=maxDifficulty)
			valData = synthDataGen.getValidation()
			if i == 0 or reset:
				saStats = sa.analyzeOneLingam(dataPoints, gen=False, valData=valData)
				difficulty = round(saStats['difficulty'],2)
				diffNorm = round(saStats['normDifficulty'],2)
				print('difficulty = ', difficulty, ', norm difficulty = ', diffNorm)
		elif testType == 'live':
			outFile = datafile
		else:
			print('*** Invalid Test Type = ', testType)
			return
		#print('Test Algorithm = ', testAlgo)
		startTime = time.time()
		method = None
		if testAlgo == 'LI':
			algo = CaLingam
		elif testAlgo == 'PL':
			algo = CaLingam
			method = 'PL'
		elif testAlgo == 'MI':
			algo = CaLingam
			method = 'MI'		
		elif testAlgo == 'PC':
			algo = CaPC
		elif testAlgo == 'CM':
			algo = CaCCM
		elif testAlgo == 'ANY':
			algo = CaICAny
		if method is not None:
			# Call algo with method parameter
			c = algo.IC(outFile, limit=dataPoints, prune=prune, method=method)
		else:
			c = algo.IC(outFile, limit=dataPoints, prune=prune)
		if i==0 and testAlgo == 'ANY':
				print('Method = ', c.method())
		dag = c.getDAG()
		endTime = time.time()
		duration = endTime - startTime
		totalDuration += duration
		if len(validation) > 0:
			result = eval(validation + '(paramSet, dag, valData)')
		else:
			corder = dag.getEquivalencyKey()
			if not corder in CAUSAL_ORDERS:
				CAUSAL_ORDERS[corder] = 1
				DAGS[corder] = dag
				#print('Causal Order = ', dag.getVarOrder())
				# If we get a new causal order after the first time, it is considered an error
				if len(CAUSAL_ORDERS) > 1:
					result = 0
				else:
					result = 1
			else:
				# Got an existing causal order.  Consider that success
				count = CAUSAL_ORDERS[corder]
				count += 1
				CAUSAL_ORDERS[corder] = count
				result = 1
		if result:
			#print ('Success', )
			print ('.',end='', flush=True)
		else:
			print ('x',end='', flush=True)
			fails += 1
	print()
	reliability = round(1.0 - (fails / float(runs)),2)
	#stats = 'Errors: ' + str(int(fails)) + ' / ' + str(int(i+1)) + ' -- Reliability: ' + str((1.0 - (fails / (i+1))) * 100) + '%'
	stats = 'Errors: ' + str(int(fails)) + ' / ' + str(int(runs)) + ' -- Reliability: ' + str(reliability * 100) + '% avg duration: ' + str(round(totalDuration/runs, 2)) + ' sec' + ' difficulty: ' + str(difficulty) + ' diffNorm: ' + str(diffNorm) + ' strength: ' + str(round(reliability * diffNorm,2))
	print('Stats = ', stats)
	if fails > 0:
		if len(validation) > 0:
			# Sort validation_errors to show the most common first
			counts = []
			keys = []
			for key in VALIDATION_ERRORS.keys():
				count = VALIDATION_ERRORS[key]
				keys.append(key)
				counts.append(count)
			tuples = list(zip(counts, keys))
			tuples.sort()
			tuples.reverse()
			keys = [key for (count, key) in tuples]
			errStrs = []
			for key in keys:
				errStrs.append(key + ':' + str(VALIDATION_ERRORS[key]))
				errStr = str.join(', ', errStrs)
			print('ValidationErrors = ', errStr)
	if len(validation) == 0:
		maxKey = None
		maxKeyCount = 0
		totalKeyCount = 0
		keys = CAUSAL_ORDERS.keys()
		for key in keys:
			keyCount = CAUSAL_ORDERS[key]
			if keyCount > maxKeyCount:
				maxKeyCount = keyCount
				maxKey = key
			totalKeyCount += keyCount
		print('Most Common Order (', maxKeyCount / totalKeyCount * 100, '%) = ', maxKey)
		print('Representative DAG = ', DAGS[maxKey])
		print()
		print('CausalOrders = ', str(CAUSAL_ORDERS))
	return


def SynthStructVal(paramSet, dag, valData):
	if valData is None:
		valData = eval(paramSet[4])
	failCount = 0
	for tuple in valData:
		successor, predecessors = tuple
		parents = dag.getParents(successor)
		for parent in parents:
			if parent not in predecessors:
				# Extra Parent
				failCount += 1
				errCount = 0
				signature = 'Extra Parent ' + parent + ' for node ' + successor + ' -- ' + str(tuple)
				if signature in VALIDATION_ERRORS:
					errCount = VALIDATION_ERRORS[signature]
				else:
					#print('Error:  ', signature)
					pass
				errCount += 1
				VALIDATION_ERRORS[signature] = errCount

		for predecessor in predecessors:
			if predecessor not in parents:
				# Missing Parent
				failCount += 1
				errCount = 0
				signature = 'Missing Parent ' + predecessor + ' for node ' + successor + ' -- ' + str(tuple)
				if signature in VALIDATION_ERRORS:
					errCount = VALIDATION_ERRORS[signature]
				else:
					#print('Error:  ', signature)
					pass
					
				errCount += 1
				VALIDATION_ERRORS[signature] = errCount
	if failCount > 0:
		return False
	return True

def SynthOrderVal(paramSet, dag, valData):
	if valData is None:
		valData = eval(paramSet[4])
	failCount = 0
	for tuple in valData:
		successor, predecessors = tuple
		ancestors = dag.getAncestors(successor)
		for predecessor in predecessors:
			if predecessor not in ancestors:
				# Missing Parent
				failCount += 1
				errCount = 0
				signature = 'Missing Ancestor ' + predecessor + ' for node ' + successor + ' -- ' + str(tuple)
				if signature in VALIDATION_ERRORS:
					errCount = VALIDATION_ERRORS[signature]
				else:
					#print('Error:  ', signature)
					pass
					
				errCount += 1
				VALIDATION_ERRORS[signature] = errCount
	if failCount > 0:
		return False
	return True
	
selections = []
if len(sys.argv) > 1:
	selectionStr = sys.argv[1]
	selTokens = selectionStr.split(',')
	for selToken in selTokens:
		selections.append(selToken.strip())
runTests(selections)
# import threading
# import time
# t = threading.Thread(target=runTests, args=(selections,))
# t.start()
# while 1:
	# time.sleep(1)
