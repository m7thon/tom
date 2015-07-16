import unittest
import tom
import os

current_dir = os.path.abspath(os.path.dirname(__file__)) + '/'

import pickle
import bz2
import random as rnd


# Provide unit type tests for the components of the tom library:
def testRandom():
    pass

def canonicalizeSuffixTree(stree):
    cForm = []
    it = tom.PrefixIterator(stree)
    while it.isValid():
        nodeData = [it.nodeIndex(), it.headIndex(), it.depth()]
        nodeData += [it.getChild().nodeIndex() if it.getChild().isValid() else -1]
        nodeData += [it.getSibling().nodeIndex() if it.getSibling().isValid() else -1]
        nodeData += [it.getSuffixLink().nodeIndex() if it.isNode() else -1]
        nodeData += [it.count()]
        cForm += [nodeData]
        it.next()
    return cForm, stree.getDeepestVirtualLeafBranch().nodeIndex()


def verifySuffixTreeCounts(seq, stree):
    it = tom.PrefixIterator(stree)
    while it.isValid():
        if seq.count(it.string()) != it.count(): return False
        it.next()
    return True


def testSuffixTree():
    print('Beginning SuffixTree tests...')
    testSeqs = ['10', 'mississippi', '10_10000', '3_2_10000', '4x3_10000', 'pc_10000']
    ok = True
    for s in testSeqs:
        print('   testing with ' + s + '.seq')
        seq = tom.load(s + '.seq')
        f = bz2.BZ2File(s + '.stree.bz2')
        streeForm = pickle.load(f)
        f.close()
        symbolSize = 2 if seq.nU() > 1 else 1
        # create SuffixTree in several random steps:
        sizes = sorted([rnd.randint(1, seq.length() - 1) for i in range(5)]) + [seq.length()]
        stree = tom.STree(seq, symbolSize, symbolSize)
        for sz in sizes:
            stree.extendTo(symbolSize * sz)
            if not verifySuffixTreeCounts(seq.sub(0, symbolSize * sz), stree):
                print('!!! SuffixTree test NOT passed !!! Substring count discrepancy.')
                ok = False
        streeFormNew = canonicalizeSuffixTree(stree)
        if not streeFormNew == streeForm:
            print('!!! SuffixTree test NOT passed !!! Something has changed.')
            ok = False
    if ok:
        print('...SuffixTree tests passed.')

# testSuffixTree()
