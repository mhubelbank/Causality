import numpy as np
exp1datafn = '..\if-dat.csv'
outfn = '..\data\experiment1.csv'
outCorrFn = '..\experiment1_cor.csv'
pthreshold = .6
ithreshold = -pthreshold

f = open(exp1datafn, 'r')

def valuesToInt(inValStrings, maxLen):
	outVals = []
	for valStr in inValStrings:
		outVals.append(int(valStr))
	while len(outVals) < maxLen:
		#print ('appending zero', len(outVals), maxLen)
		outVals.append(0)
	if len(outVals) > maxLen:
		print ('more than maxLen values ',len(outVals), maxLen)
	elif len(outVals) < maxLen:
		print ('less than maxLen values ', len(outVals), maxLen)
	return outVals

NameMap = {'138.42.94.173|Se2/0/3|bytesIn':'Ai',
			'138.42.94.173|Se2/0/3|bytesOut':'Ao',
			'138.42.94.173|Se2/0/1|bytesIn':'Bi',
			'138.42.94.173|Se2/0/1|bytesOut':'Bo',
			'138.42.96.6|Fa1/0|bytesIn':'Ci',
			'138.42.96.6|Fa1/0|bytesOut':'Co',
			'138.42.96.5|Fa3/0|bytesIn':'Hi',
			'138.42.96.5|Fa3/0|bytesOut':'Ho',
			'138.42.94.197|Et0/5|bytesIn':'Gi',
			'138.42.94.197|Et0/5|bytesOut':'Go',
			'138.42.94.173|Fa3/0/0|bytesIn':'Di',
			'138.42.94.173|Fa3/0/0|bytesOut':'Do',
			'138.42.94.173|Fa3/1/0|bytesIn':'Ei',
			'138.42.94.173|Fa3/1/0|bytesOut':'Eo',
			'138.42.96.6|Fa3/0|bytesIn':'Fi',
			'138.42.96.6|Fa3/0|bytesOut':'Fo',
			'138.42.96.5|Gi2/0|bytesIn':'Ii',
			'138.42.96.5|Gi2/0|bytesOut':'Io',
			'138.42.96.5|Fa1/0|bytesIn':'Ji',
			'138.42.96.5|Fa1/0|bytesOut':'Jo',
			}

def run():
	maxsamples = 0
	devlist = []
	#devlist = ['138.42.94.197', '138.42.96.5', '138.42.96.6']
	#devlist = ['138.42.96.6','138.42.96.36', '138.42.94.173']
	devlist = ['138.42.96.6', '138.42.94.173', '138.42.96.5', '138.42.94.197']
	
	#iflist = ['Fa1/0','Gi0/0/0','Gi1/0/1','Fa3/0/0','Gi3/25','Se2/0/1','Fa4/0','Fa2/2/9']
	iflist = ['Fa1/0', 'Fa3/0', 'Fa3/0/0', 'Se2/0/3', 'Se2/0/1', 'Et0/5','Fa3/1/0', 'Gi2/0']
	dupdict = {}
	d = {}
	print ('f = ', f)
	count = 0
	sorted = []
	while 1:
		l = f.readline()
		if len(l) <= 0:
			break
		
		sorted.append(l)
	sorted.sort()
	lastSampDict = {}
	for l in sorted:
		if count > 100000000:
			break
		tokens = l.split(',')
		#print ('t=', tokens)
		
		if len(tokens) >= 8:
			timestamp, time, devname, ifname, bytesIn, bytesOut, errorsIn, errorsOut = tokens[:8]
			duptag = timestamp + '|' + devname + '|' + ifname
			if duptag in dupdict:
				continue
			if len(devlist) > 0 and devname not in devlist:
				continue
			if len(iflist) > 0 and ifname not in iflist:
				continue
			dupdict[duptag] = 1
			tag = devname + '|' + ifname + '|bytesIn'
			if tag not in d:
				d[tag] = (tag, ['0']*(maxsamples-1))
			dtag, vals = d[tag]
			vlen = len(vals)
			while vlen < maxsamples - 1:
				vals.append(lastSampDict[tag])
				vlen += 1
			lastSampDict[tag] = bytesIn
			vals.append(bytesIn)
			sampCount = len(vals)
			if sampCount > maxsamples:
				maxsamples = sampCount
				#print ('maxsamples1 =', timestamp, maxsamples, count, len(d), tag, vals)
			tag = devname + '|' + ifname + '|bytesOut'
			if tag not in d:
				d[tag] = (tag, ['0']*(maxsamples-1))
			dtag, vals = d[tag]
			vlen = len(vals)
			while vlen < (maxsamples - 1):
				vals.append(lastSampDict[tag])
				vlen += 1
			lastSampDict[tag] = bytesOut
			vals.append(bytesOut)
			sampCount = len(vals)
			if sampCount > maxsamples:
				maxsamples = sampCount
				#print ('maxsamples2 =', timestamp, maxsamples, count, len(d), tag, value)
		count += 1
	print ('rawcount=', count)
	
	arrayData = []
	arrayTags = []
	keys = list(d.keys())
	keys.sort()
	outlines = []
	count = 0
	sampcount = 1000000
	sampstart = 0
	for key in keys:
		if count==500:
			break
		entry = d[key]
		values = entry[1]
		nonzero = False
		tokens = str.split(key, '|')
		#print ('tokens = ', str(tokens))
		ifname = tokens[1]
		iftokens = str.split(ifname, '/')
		if len(iftokens) > 1:
			pass
			#continue
		iftokens = str.split(ifname, '.')
		if len(iftokens) > 1:
			pass
			#continue
		for value in values:
			if value != '0':
				nonzero = True
		if nonzero:
			#outline = entry[0] + ',' + str.join(',',values) + '\n'
			#outlines.append(outline)
			arrayTags.append(key)
			arrayData.append(valuesToInt(values, maxsamples))
			#print ('valueslen = ', len(values), len(valuesToInt(values, maxsamples)))
			count += 1
	#print ('arraydata = ', arrayData)
	a = np.array(arrayData)
	samples = min(maxsamples, sampcount)
	f2 = open(outfn, 'w')
	for i in range(sampstart, sampstart+samples):
		row = []
		for j in range(len(arrayData)):
			row.append(str(arrayData[j][i]))
		outline = str.join(',', row)
		if i < sampstart+samples - 1:
			outline += '\n'
		outlines.append(outline)
	transTags = []
	for tag in arrayTags:
		try:
			transTag = NameMap[tag]
		except:
			print('Tag Translation Failed: ', tag)
			transTag = tag
		transTags.append(transTag)
	print('Raw Tags are: ', arrayTags)
	print('Translated Tags are: ', transTags)
	f2.writelines([str.join(',', transTags)+'\n'])
	f2.writelines(outlines)
	f2.close()

	print ('a=', a.shape)
	cova = np.corrcoef(a)
	
	#print ('cov  = ', cova)
	isCorrelated = []
	for val in cova.flat:
		if val > pthreshold or val < ithreshold:
			isCorrelated.append(1)
		else:
			isCorrelated.append(0)
	corrmat = np.zeros(cova.shape, np.int8)
	corrmat.put(np.arange(len(isCorrelated)), isCorrelated)
	corrdict = {}
	corrlen, corrwidth = corrmat.shape
	finalCorrCount = 0
	throttlecount = 0
	if throttlecount > 0:
		if corrlen > throttlecount:
			corrlen = throttlecount
			corrwidth = throttlecount
	for i in range(corrlen):
		name1 = arrayTags[i]
		for j in range(corrwidth):
			if i == j:
				continue
			if corrmat[i,j] == 1:
				finalCorrCount += 1
				name2 = arrayTags[j]
				val = cova[i][j]
				print (name1 + '(' + str(i) + ')', '<==>', name2 + '(' + str(j) + ') -- ', val)
	print ('cellCount = ', corrlen * corrwidth / 2)
	print ('finalCorrCount = ', finalCorrCount)
	print ('corrPercent = ', (finalCorrCount*100) /(corrlen * corrwidth / 2))
	#print ('isCorrelated', isCorrelated)
	#print ('corrmat=', corrmat)
	#print ('outlines = ', len(outlines))
	# Generate correlation matrix file
	corrFile = open(outCorrFn, 'w')
	lines = []
	# Header
	linetokens = []
	for i in range (corrlen):
		linetokens.append(arrayTags[i])
	corrFile.write(str(samples) + '\n')
	corrFile.write(str.join(',', linetokens) + '\n')
	# Now the body
	for i in range (corrlen):
		linetokens = []
		for j in range(corrwidth):
			if j > i:
				break
			tag1 = arrayTags[i]
			tag2 = arrayTags[j]
			tokens1 = tag1.split('|')
			tokens2 = tag2.split('|')
			el1 = str.join('|', tokens1[:-1])
			el2 = str.join('|', tokens2[:-1])
			v1 = tokens1[-1]
			v2 = tokens2[-1]
			coef = str(cova[i,j])
			if (v1 == v2) and (el1 != el2): # No in to in or out to out
				coef = '0'
			elif (el1 == el2) and (v1 != v2): # No in / out for same element
				coef = '0'
			linetokens.append(coef)
		if i < corrlen-1:
			corrFile.write(str.join(',', linetokens) + '\n')
		else:
			corrFile.write(str.join(',', linetokens))
	corrFile.close()
	
			
run()
