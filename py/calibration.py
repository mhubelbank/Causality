

class Calibration:
	def __init__(self, clsName):
		self.clsName  = clsName
		self.vars = {}
		filename = clsName + '.cal'
		try:
			f = open(filename, 'r')
			instr = f.read()
			f.close()
		except:
			instr = ''
		lines = instr.split('\n')
		for line in lines:
			if len(line) == 0:
				continue
			tokens = line.split('=')
			if len(tokens) != 2:
				continue
			varName = tokens[0].strip()
			varVal = tokens[1].strip()
			self.vars[varName] = varVal
		return
		
	def get(self, vName, default=None):
		return self.vars.get(vName, default)
		
	def set(self, vName, vValue):
		self.vars[vName] = vValue
		return
		
	def save(self,*args):
		filename = self.clsName + '.cal'
		f = open(filename, 'w')
		outLines = []
		for vName in self.vars.keys():
			vValue = self.vars[vName]
			outLine = vName + ' = ' + vValue
			outLines.append(outLine)
		outStr = str.join('\n', outLines)
		f.write(outStr)
		f.close()
	
		