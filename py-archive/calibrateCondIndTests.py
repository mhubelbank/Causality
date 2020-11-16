import analyzeSEM
import CondIndGauss
import CondIndGCR
import getData
import synthDataGen
import calibration
import sys


FILES = ['..\Tests\synthTest1', '..\Tests\synthTest2', '..\Tests\synthTest2a', '..\Tests\synthTest4', '..\Tests\synthTest6', '..\Tests\synthTest8']
FILES = ['..\Tests\synthTest1', '..\Tests\synthTest2', '..\Tests\synthTest2a', '..\Tests\synthTest8']
#FILES = ['..\Tests\synthTest2']
TEST_TYPES = ['CondIndGauss', 'CondIndGCR']
TEST_TYPES = ['CondIndGCR']
ERROR_FREE_RUNS = 10 # OF error-free runs needed to declare victory

class Calibrator:
	def __init__(self, runs=100, datacount=10000, noSave=False):
		self.datacount = datacount
		for testType in TEST_TYPES:
			print('Calibrating Independence Test: ', testType)
			cal = calibration.Calibration(testType)
			self.origThreshold = eval(testType + '.DEPENDENCE_THRESHOLD')
			self.currThreshold = self.origThreshold
			prevDirection = 0
			#scale = .1
			flipFlops = 2
			zeroCount = 0
			for i in range(runs):
				self.testCount = 0
				self.err1Count = 0 # Type 1 Error.  Should be dependent, but was found independent
				self.err2Count = 0 # Type 2 Error.  Should be independent, but was found dependent
				for filePath in FILES:
					self.calibrateOneCITest(testType, filePath)
				errors = self.err1Count + self.err2Count
				direction = self.err2Count - self.err1Count
				print('Score = ', (1 - (errors / self.testCount))*100, '%   ', self.err1Count, self.err2Count, direction)
				if errors == 0:
					print('No Errors -- Perfect Calibration')
					zeroCount += 1
					if zeroCount > ERROR_FREE_RUNS:
						break
				else:
					zeroCount = 0
				if flipFlops > 16:
					break
				if (prevDirection > 0 and direction < 0) or (prevDirection < 0 and direction > 0):
					flipFlops += 1
				#print('flipFlops = ', flipFlops, prevDirection, direction)
				if direction < 0:
					adjustment = -self.origThreshold * 2**-(flipFlops)
				elif direction > 0:
					adjustment = self.origThreshold * 2**-(flipFlops)
				else:
					print('Errors are balanced')
					adjustment = 0.0
				if adjustment != 0:
					self.currThreshold += adjustment
					print('Adjusting threshold by ', adjustment, 'to', self.currThreshold)
					exec(testType + '.DEPENDENCE_THRESHOLD = ' + str(self.currThreshold))
					cal.set('DEPENDENCE_THRESHOLD', str(self.currThreshold))
					#cal.save()
					prevDirection = direction
			if not noSave:
				cal.save()
			print('Finished Calibration for ', testType)
		return
		
	def calibrateOneCITest(self, testType, filePath):
		synthDataGen.run(filePath + '.py', samples=self.datacount)
		exec('import ' + testType)
		module = eval(testType)
		SA = analyzeSEM.SemAnalyzer(filePath + '.py', self.datacount)
		reader = getData.DataReader(filePath + '.csv', self.datacount)
		dependencies, independencies = SA.getCondDependencies()
		# print('dependencies = ', dependencies)
		# print('independencies = ', independencies)
		errors = 0
		errorTerms={}
		items = 0
		for item in dependencies:
			x, y, z = item
			X = reader.getSeries(x)
			Y = reader.getSeries(y)
			Z = reader.getSeries(z)
			ind = module.isIndependent(X,Y, Z)
			if ind:
				print('Error -- ', x, 'and', y, 'Should be dependent given', z)
				self.err1Count += 1
				errors += 1
				errorTerms[item]=1
			self.testCount += 1
		for item in independencies:
			x, y, z = item
			X = reader.getSeries(x)
			Y = reader.getSeries(y)
			Z = reader.getSeries(z)
			ind = module.isIndependent(X,Y,Z)
			if not ind:
				print('Error -- ', x, 'and', y, 'Should be independent given', z)
				self.err2Count += 1
				errors += 1
				errorTerms[item] = 1
			self.testCount += 1
		#print('Rating = ', (1 - (errors / items))*100, '%')
		print('Errors for file: ',filePath, '=', errors, list(errorTerms.keys()))
		return

noSave = True		
if 'save' in sys.argv:
	noSave = False
	print('Save is enabled')
else:
	print('Save is disabled')
c = Calibrator(datacount=1000, noSave=noSave)
