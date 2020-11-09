import sys

DEBUG=False
if '-D' in sys.argv:
	sys.argv.remove('-D')
	DEBUG=True

def Dprint(*args):
	if DEBUG:
		print(*args)

