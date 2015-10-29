from __future__ import (division, absolute_import, print_function, unicode_literals)
from .. import _tomlib

def stringToSequence(str, symbolTable=None):
    if not symbolTable: symbolTable = {}
    char = set()
    for c in str:
        if c not in symbolTable:
            char.add(c)
    char = sorted(list(char), key=lambda c: str.count(c), reverse=True)
    firstFreeSymbolValue = 0 if len(symbolTable) == 0 else max(symbolTable.values()) + 1
    for i in range(len(char)):
        symbolTable[char[i]] = firstFreeSymbolValue + i
    return _tomlib.Sequence([symbolTable[c] for c in str], firstFreeSymbolValue + len(char)), symbolTable
