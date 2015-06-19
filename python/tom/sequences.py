from .tomlib import *

def stringToSequence(str, symbolTable=None):
	if not symbolTable: symbolTable = {}
	char = set()
	for c in str:
		if c not in symbolTable:
			char.add(c)
	char = sorted(list(char), key= lambda c: str.count(c), reverse=True)
	firstFreeSymbolValue = 0 if len(symbolTable) == 0 else max(symbolTable.values()) + 1
	for i in range(len(char)):
		symbolTable[char[i]] = firstFreeSymbolValue + i
	return Sequence([symbolTable[c] for c in str], firstFreeSymbolValue + len(char)), symbolTable

def incrementSequenceList(sl, nO, nU=0):
    io = (nU != 0)
    nUO = (nU, nO)
    index = len(sl) - 1
    carry = True
    while carry:
        if index < 0:
            if io:
                sl.append(0)
            sl.append(1)
            return sl
        sl[index] += 1
        if (io and sl[index] == nUO[index%2]) or (not io and sl[index] == nO):
            sl[index] = 0
            index -= 1
        else:
            carry = False
    return sl

def generateSequences(len, nO, nU=0):
    """Generate the set of all Sequences of length len in lexicographic order
    with input alphabet size nU (which defaults to 0) and output alphabet size nO."""
    io = (nU != 0)
    s = len * [0] if not io else len * [0,0]
    seqs = Sequences()
    for i in range((nO * max(nU,1))**len):
        seqs.push_back(Sequence(s, nO, nU))
        s = incrementSequenceList(s, nO, nU)
    return seqs
