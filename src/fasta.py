import re

from .sequence import Sequence

def readFasta(file):
	with open(file, 'r') as f:
		content = ''.join(f.readlines())
		m = re.search('>(.*)\n([\s\S]+?(?=>|$))', content) # find the first fasta input
		return m[0][1:], Sequence(m[1], m[2].replace('\n', ''))

def readFastaMul(file):
	with open(file, 'r') as f:
		content = ''.join(f.readlines())
		return tuple(Sequence(m.group(1), m.group(2).replace('\n','')) for m in re.finditer('>(.*)\n([\s\S]+?(?=>|$))', content))
