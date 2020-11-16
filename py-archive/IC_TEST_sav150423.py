# IC_TEST
#
# Master Test for the Inferred Causality (IC) Prototype Subsystem

# PARAMS are: (<datafile>, <runs>, <validationFunction>)
from debug import *
import CaLingam
import synthDataGen
import sys
import time
import analyzeSEM

VALIDATION_ERRORS = {}
CAUSAL_ORDERS = {}

TESTS = [
		# TestDef := (<testName>, <testType := 'synth' | 'live'>, <testDescription>, <regenCount or 0>, <maxDifficulty or 0>, [<Params>])
		# regenCount :=  Number of times to run before resetting the data generator to create a new parameterization.  Zero to disable.
		# maxDifficulty := The maximum difficulty for the data generator to generate or zero to disable.  Note: This can cause a reduction in
		#    performance of the tests, or even stalling if attempting to create tests with too small a maxDifficulty.

		('Test1', 'synth', 'Basic A->I Pattern', 0, 0, [
			# Params := (<DataGeneratorFile(synth) or DataFile(live)>, <how many times to run>, <how many data points per variable>, <validation function>, <validation data override(optional-->uses validation from input file by default)

			('..\Tests\synthTest1.py', 10, 10000, 'SynthOrderVal'),
			('..\Tests\synthTest1.py', 100, 500, 'SynthStructVal'),
			('..\Tests\synthTest1.py', 100, 100, 'SynthStructVal'),
			('..\Tests\synthTest1.py', 100, 50, 'SynthStructVal'),
			('..\Tests\synthTest1.py', 100, 10, 'SynthStructVal'),
			]),
		('Test2', 'synth', 'Extended A->L Pattern', 0, 0, [
			('..\Tests\synthTest2.py', 100, 1000, 'SynthOrderVal'),
			('..\Tests\synthTest2.py', 100, 500, 'SynthStructVal'),
			('..\Tests\synthTest2.py', 100, 100, 'SynthStructVal'),
			('..\Tests\synthTest2.py', 100, 50, 'SynthStructVal'),
			('..\Tests\synthTest2.py', 100, 10, 'SynthStructVal'),
			]),
		('Test3', 'synth', 'Basic A->I with non-linear terms', 0, 0, [
			('..\Tests\synthTest3.py', 100, 1000, 'SynthOrderVal'),
			('..\Tests\synthTest3.py', 100, 500, 'SynthStructVal'),
			('..\Tests\synthTest3.py', 100, 100, 'SynthStructVal'),
			('..\Tests\synthTest3.py', 100, 50, 'SynthStructVal'),
			('..\Tests\synthTest3.py', 100, 10, 'SynthStructVal'),
			]),
		('Test4', 'synth', 'Basic A->I with simulated time-series i.e., non-iid', 0, 0, [
			('..\Tests\synthTest4.py', 100, 10000, 'SynthOrderVal'),
			('..\Tests\synthTest4.py', 100, 500, 'SynthStructVal'),
			('..\Tests\synthTest4.py', 100, 100, 'SynthStructVal'),
			('..\Tests\synthTest4.py', 100, 50, 'SynthStructVal'),
			('..\Tests\synthTest4.py', 100, 10, 'SynthStructVal'),
			]),
		('Test5', 'synth', 'Experimental', 0, 0, [
			('..\Tests\synthTest5.py', 20, 1000, 'SynthStructVal'),
			]),
			
		('Test6', 'synth', 'A->I with random parameterization', 10, 2.5, [
			('..\Tests\synthTest6.py', 100, 1000, 'SynthOrderVal'),
			('..\Tests\synthTest6.py', 100, 8000, 'SynthOrderVal'),
			('..\Tests\synthTest6.py', 100, 5000, 'SynthOrderVal'),
			('..\Tests\synthTest6.py', 100, 3000, 'SynthOrderVal'),
			('..\Tests\synthTest6.py', 100, 1000, 'SynthOrderVal'),		
			]),
		('Test7', 'synth', 'Two variable test with random parameterization', 10, 25.0, [
			('..\Tests\synthTest7.py', 100, 100, 'SynthOrderVal'),
			('..\Tests\synthTest7.py', 100, 200, 'SynthOrderVal'),
			('..\Tests\synthTest7.py', 100, 500, 'SynthOrderVal'),
			('..\Tests\synthTest7.py', 100, 1000, 'SynthOrderVal'),
			('..\Tests\synthTest7.py', 100, 2000, 'SynthOrderVal'),
			('..\Tests\synthTest7.py', 100, 5000, 'SynthOrderVal'),
			]),
		('Test8', 'synth', 'Dougs A->I test with all coefs and stddevs = 0', 0, 0, [
			('..\Tests\synthTest8.py', 10, 10000, 'SynthOrderVal'),
			('..\Tests\synthTest8.py', 10, 8000, 'SynthOrderVal'),
			('..\Tests\synthTest8.py', 10, 5000, 'SynthOrderVal'),
			('..\Tests\synthTest8.py', 10, 3000, 'SynthOrderVal'),
			('..\Tests\synthTest8.py', 10, 1000, 'SynthOrderVal'),
			]),
		('Test9', 'synth', 'Experimental', 0, 0, [
			('..\Tests\synthTest9.py', 100, 20000, 'SynthOrderVal'),
			('..\Tests\synthTest9.py', 200, 10000, 'SynthOrderVal'),
			('..\Tests\synthTest9.py', 400, 5000, 'SynthOrderVal'),
			('..\Tests\synthTest9.py', 600, 2000, 'SynthOrderVal'),
			('..\Tests\synthTest9.py', 800, 1000, 'SynthOrderVal'),
			('..\Tests\synthTest9.py', 1000, 100, 'SynthOrderVal'),
			]),
		('Test10', 'synth', 'Non-alphabetic order test', 0, 0, [
			('..\Tests\synthTest10.py', 10, 10000, 'SynthOrderVal'),
			]),
		('OneRouterTest', 'live', 'Live data with single router interface', 0, 0, [
			('..\experiment1.csv', 100, 1000, ''),
			]),
		]
		
def runTests(selectionList):
	global RESET_COUNT, MAX_DIFFICULTY
	for test in TESTS:
		testName, testType, testDescr, resetCount, maxDifficulty, testParams = test[:6]
		RESET_COUNT = resetCount
		MAX_DIFFICULTY = maxDifficulty
		if selectionList is None or len(selectionList) == 0 or testName in selectionList:
			print('Starting Test: ', testName, '(', testDescr, ')')
			runTest(testName, testType, testDescr, testParams)
			print('Test: ', testName, ' Complete\n')
	return
	
def runTest(testName, testType, testDescr, testParams):
		
	for paramSet in testParams:
		# Run through one test
		#print()
		print('Running test with parms = ', paramSet)
		datafile, runs, validation = paramSet[:3]
		CaLingamTest(testType, paramSet)
	return

	
RESET_COUNT = 0
MAX_DIFFICULTY = 0

def CaLingamTest(testType, paramSet):
	global VALIDATION_ERRORS, CAUSAL_ORDERS
	VALIDATION_ERRORS = {}
	CAUSAL_ORDERS = {}
	maxDifficulty = MAX_DIFFICULTY
	fails = 0.0
	datafile, runs, dataPoints, validation = paramSet[:4]
	valData = None
	totalDuration = 0
	datafileRootName = datafile.split('.')[-2].split('\\')[-1]
	#print('dfrn = ', datafileRootName)
	sa = analyzeSEM.SemAnalyzer(datafileRootName)
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
				saStats = sa.analyzeOne(dataPoints, gen=False, valData=valData)
				difficulty = round(saStats['difficulty'],2)
				diffNorm = round(saStats['normDifficulty'],2)
				print('difficulty = ', difficulty, ', norm difficulty = ', diffNorm)
		elif testType == 'live':
			outFile = datafile
		else:
			print('*** Invalid Test Type = ', testType)
			return
		startTime = time.time()
		c = CaLingam.IC(outFile, limit=dataPoints, prune=prune)
		dag = c.getDAG()
		endTime = time.time()
		duration = endTime - startTime
		totalDuration += duration
		if len(validation) > 0:
			result = eval(validation + '(paramSet, dag, valData)')
		else:
			corder = str.join('||', dag.getVarOrder()[:5])
			if not corder in CAUSAL_ORDERS:
				CAUSAL_ORDERS[corder] = 1
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
		else:
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
			print()
			print('CausalOrders = ', str(CAUSAL_ORDERS))
	return

calValData1 = [('E', ['A','B','C','D']),('A', []), ('B', []), ('F', ['E']), ('G', ['F']), ('H',['F']),	('I', ['F']), 
					]
calValData2 = [('E', 'A'),('E', 'B'), ('E', 'C'), ('E','D'), ('F', 'E'), ('G', 'F'), ('H','F'), 
					('I', 'F'), ('K', 'C'), ('K', 'D'), ('L', 'K'), ('L','I'), ('L','J'),
					]
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
