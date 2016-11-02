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

def verifySuffixTreeForm(seq, cform1, cform2):
    if len(cform1) != len(cform2): return False
    if cform1[-1] != cform2[-1]: return False
    for n1, n2 in zip(cform1[0], cform2[0]):
        if n1[0] != n2[0] or n1[2:] != n2[2:]: return False
        if n1[1] != n2[1] and seq.rawSub(n1[1],n1[2]) != seq.rawSub(n2[1],n2[2]): return False
    return True

class TestSuffixTree(unittest.TestCase):
    def test_STree(self):
        testSeqs = ['10', 'mississippi', '10_10000', '4_3_10000', '4x3_10000', 'pc_10000']
        for s in testSeqs:
            print('.', end='', flush=True)
            seq = tom.load(current_dir + s + '.seq.bz2')
            f = bz2.BZ2File(current_dir + s + '.stree.bz2')
            streeForm = pickle.load(f)
            f.close()
            symbolSize = 2 if seq.nInputSymbols() > 0 else 1
            # create SuffixTree in several random steps:
            sizes = sorted([rnd.randint(1, seq.length() - 1) for i in range(5)]) + [seq.length()]
            stree = tom.STree(seq.sub(0,0))
            for sz in sizes:
                stree.extendTo(seq.sub(0,sz))
                self.assertTrue(verifySuffixTreeCounts(seq.sub(0, sz), stree),
                                "!!! SuffixTree test NOT passed !!! Substring count discrepancy.")
            streeFormNew = canonicalizeSuffixTree(stree)
            self.assertTrue(verifySuffixTreeForm(seq, streeFormNew, streeForm), "!!! SuffixTree test NOT passed !!! Something has changed.")
            for l in range(1, min(seq.length(), 12)):
                self.assertTrue(tom.stree.Position(stree, seq.slice(-l)).isSuffix(), ".isSuffix() test failed!")

        print(' ', end='', flush=True)
