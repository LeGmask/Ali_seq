import re
from typing import Tuple

from .sequence import Sequence

def readFasta(file: str) : #-> Tuple(str, Sequence)
	with open(file, 'r') as f:
		content = ''.join(f.readlines())
		m = re.search('>(.*)\n([\s\S]+?(?=>|$))', content) # find the first fasta input
		return m[0][1:], Sequence(m[1], m[2].replace('\n', ''))

def readFastaMul(file: str) : # -> Tuple(Sequence)
	with open(file, 'r') as f:
		content = ''.join(f.readlines())
		return tuple(Sequence(m.group(1), m.group(2).replace('\n','')) for m in re.finditer('>(.*)\n([\s\S]+?(?=>|$))', content))
