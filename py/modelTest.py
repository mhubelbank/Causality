# model test
import CaLingam
import synthDataGen
import sys

# Parameters for getBestModel
NUMBER_OF_MODELS_TO_CONSIDER = 11 # Pick most common model out of N tries

def run(filePath, dataPoints):
	model = getBestModel(filePath, dataPoints)
	model.calculate()
	print('Original model: ', str(model))
	# xNodes = model.getExogenousNodes()
	# print('xNodes = ', xNodes)
	# testNode = xNodes[0]
	# origValue = model.getNodeValue(testNode)
	# origMean, origStd = origValue
	# newMean = origMean + 5
	# newStd = origStd
	# print('Adding 5 to original mean ', origMean, ' for node ', testNode)
	# model.setNodeValue(testNode, (newMean, newStd))
	# model.calculate()
	# print ('model = ', str(model))
	print('\n\nAnalysis of this model:\n')
	print(model.implications())
	return

def getBestModel(filePath, dataPoints):
	dot_pos = filePath.rfind('.')
	basePath = filePath[:dot_pos]
	ext = filePath[dot_pos+1:]
	modelDict = {}
	for i in range(NUMBER_OF_MODELS_TO_CONSIDER):
		if ext == 'py':
			#print('Generating data')
			synthDataGen.run(filePath, samples=dataPoints)
		dataFile = basePath + '.csv'
		
		ic = CaLingam.IC(dataFile, dataPoints)
		model = ic.getModel()
		dagStr = model.getEquivalencyKey()
		if not dagStr in modelDict:
			modelDict[dagStr] = []
		modelDict[dagStr].append(model)
	# Find model based on most common structure
	mostCommon = None
	mostCommonCount = 0
	for modelList in modelDict.values():
		if len(modelList) > mostCommonCount:
			mostCommonCount = len(modelList)
			mostCommon = modelList[0]
	model = mostCommon
	print('Most common model = ', mostCommonCount, '/', NUMBER_OF_MODELS_TO_CONSIDER)
	return model
	
if len(sys.argv) > 1:
	inputFile = sys.argv[1]
	if len(sys.argv) > 2:
		dataPoints = int(sys.argv[2])
	else:
		dataPoints = 1000
	run(inputFile, dataPoints)
