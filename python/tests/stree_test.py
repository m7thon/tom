import unittest
import tom
import os

current_dir = os.path.abspath(os.path.dirname(__file__)) + '/'

import pickle
import bz2
import random as rnd

def canonicalizeSuffixTree(stree):
    cForm = []
    it = tom.tomlib.PrefixIterator(stree)
    while it.isValid():
        nodeData = [it.nodeIndex(), it.headIndex(), it.depth()]
        nodeData += [it.child().nodeIndex() if it.child().isValid() else -1]
        nodeData += [it.sibling().nodeIndex() if it.sibling().isValid() else -1]
        nodeData += [it.suffixlink().nodeIndex() if it.isInternal() else -1]
        nodeData += [it.count()]
        cForm += [nodeData]
        it.next()
    return cForm, stree.getDeepestVirtualLeafBranch().nodeIndex()


def verifySuffixTreeCounts(seq, stree):
    it = tom.tomlib.PrefixIterator(stree)
    while it.isValid():
        if seq.count(it.asSequence()) != it.count():
            return False
        it.next()
    return True

class TestSuffixTree(unittest.TestCase):
    def test_STree(self):
        testSeqs = ['10', 'mississippi', '10_10000', '4_3_10000', '4x3_10000', 'pc_10000']
        print()
        for s in testSeqs:
            print('   testing with ' + s + '.seq.bz2')
            seq = tom.load(current_dir + s + '.seq.bz2')
            f = bz2.BZ2File(current_dir + s + '.stree.bz2')
            streeForm = pickle.load(f)
            f.close()
            symbolSize = 2 if seq.nU() > 0 else 1
            # create SuffixTree in several random steps:
            sizes = sorted([rnd.randint(1, seq.length() - 1) for i in range(5)]) + [seq.length()]
            stree = tom.tomlib.STree(seq, 1)
            for sz in sizes:
                stree.extendTo(sz)
                self.assertTrue(verifySuffixTreeCounts(seq.sub(0, sz), stree),
                                "!!! SuffixTree test NOT passed !!! Substring count discrepancy.")
            streeFormNew = canonicalizeSuffixTree(stree)
            self.assertTrue(streeFormNew == streeForm, "!!! SuffixTree test NOT passed !!! Something has changed.")
