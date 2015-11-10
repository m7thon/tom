import unittest
import tom
import os

current_dir = os.path.abspath(os.path.dirname(__file__)) + '/'

import pickle
import bz2
import random as rnd

def nodeIndex(node):
    return node.nidx() & (tom.stree.INTERNAL | tom.stree.INDEX)

def canonicalizeSuffixTree(stree):
    cForm = []
    for it in tom.stree.PrefixIterator(stree):
        nodeData = [nodeIndex(it), it.headIndex(), it.depth()]
        nodeData += [nodeIndex(it.child()) if it.child().isValid() else -1]
        nodeData += [nodeIndex(it.sibling()) if it.sibling().isValid() else -1]
        nodeData += [nodeIndex(it.suffix()) if it.isInternal() else -1]
        nodeData += [it.count()]
        cForm += [nodeData]
    return cForm, nodeIndex(tom.stree.Node(stree, stree.deepestInternalSuffixNidx()))


def verifySuffixTreeCounts(seq, stree):
    for node in tom.stree.PrefixIterator(stree):
        if seq.count(node.sequence()) != node.count():
            return False
    return True

class TestSuffixTree(unittest.TestCase):
    def test_STree(self):
        testSeqs = ['10', 'mississippi', '10_10000', '4_3_10000', '4x3_10000', 'pc_10000']
        print()
        for s in testSeqs:
            print('   testing with ' + s + '.seq.bz2')
            seq = tom.io.load(current_dir + s + '.seq.bz2')
            f = bz2.BZ2File(current_dir + s + '.stree.bz2')
            streeForm = pickle.load(f)
            f.close()
            symbolSize = 2 if seq.nInputSymbols() > 0 else 1
            # create SuffixTree in several random steps:
            sizes = sorted([rnd.randint(1, seq.length() - 1) for i in range(5)]) + [seq.length()]
            stree = tom.stree.STree(seq.sub(0,0))
            for sz in sizes:
                stree.extendTo(seq.sub(0,sz))
                self.assertTrue(verifySuffixTreeCounts(seq.sub(0, sz), stree),
                                "!!! SuffixTree test NOT passed !!! Substring count discrepancy.")
            streeFormNew = canonicalizeSuffixTree(stree)
            self.assertTrue(streeFormNew == streeForm, "!!! SuffixTree test NOT passed !!! Something has changed.")
