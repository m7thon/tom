import unittest
import tom
import sys,os

dir = os.path.abspath(os.path.dirname(__file__)) + '/'

def SequenceEqualsList(tom_seq, seq):
  if tom_seq.size() != len(seq): return False
  for i in range(len(seq)):
    if tom_seq[i] != seq[i]: return False
  return True

class TestSequence(unittest.TestCase):
  def test_slicing(self):
    for seq in [[], [1], [1,2], [1,2,3], list(range(7)), list(range(8))]:
      for reverse in [False, True]:
        if reverse: seq.reverse()
        for tom_seq in [tom.Sequence(seq, 10, 0), tom.Sequence(seq, 10, 1)]:
          if reverse: tom_seq.reverse()
          for i in list(range(-15,15)):
            for j in list(range(-15,15)):
              for d in [-1,1]:
                seq_slice = seq[None if i == -15 else i:None if j == -15 else j:d]
                tom_seq_slice = tom_seq[None if i == -15 else i:None if j == -15 else j:d]
                tom_slice = tom_seq.slice(tom.NoIndex if i == -15 else i, tom.NoIndex if j == -15 else j, d == -1)
                self.assertTrue(tom_seq_slice.nU() == tom_slice.nU() == tom_seq.nU(), "slicing changed the input alphabet for [%d:%d:%d]" %(i,j,d))
                self.assertTrue(tom_seq_slice.nO() == tom_slice.nO() == tom_seq.nO(), "slicing changed the output alphabet for [%d:%d:%d]" %(i,j,d))
                self.assertTrue(SequenceEqualsList(tom_slice, seq_slice), "python slicing error for [%d:%d:%d]: " %(i,j,d) + str(tom_slice) + " != " + str(seq_slice) + " case:" + str(seq))
                self.assertTrue(SequenceEqualsList(tom_seq_slice, seq_slice), "slice() error for [%d:%d:%d]" %(i,j,d) + str(tom_seq_slice) + " != " + str(seq_slice) + " case:" + str(seq))
                self.assertTrue(tom_seq_slice.copy() == tom_seq_slice, "slicing and copying does not play together for [%d:%d:%d]" %(i,j,d) + " case:" + str(tom_seq_slice))

                try:
                  tom_seq[None if i == -15 else i:None if j == -15 else j:d] = seq[None if i == -15 else i:None if j == -15 else j:d]
                  self.assertTrue(False, "slicing must not be useable as assignment")
                except TypeError:
                  pass

  def test_subindexinx(self):
    pass

  def test_equality(self):
    for l in range(5):
      for o in range(2):
        for u in range(2):
          try:
            s = tom.Sequence(l, o, u)
            self.assertTrue(not (l > 0 and o == 0), "don't allow a Sequence to be created with nO=0 but size>0")
          except:
            pass
          self.assertTrue(bool(s) == (s.nO() != 0) , "A sequence must be true iff nO != 0")
          for seq in [tom.Sequence(0,1,0), tom.Sequence(0,1,3), tom.Sequence(2, 4, 0), tom.Sequence(2, 0, 0), tom.Sequence(3, 4, 4)]:
            self.assertTrue((s == seq) == (s.size() == seq.size() and s.nO() != 0 and seq.nO() != 0 and
                                           (s.nU() != 0) == (seq.nU() != 0)), str(s) + "==" + str(seq))

  def test_json_io(self):
    for seq in [[], [1], [1,2], [1,2,3], list(range(7)), list(range(8))]:
      for tom_seq in [tom.Sequence(seq, 10, 0), tom.Sequence(seq, 10, 1)]:
        dummy()


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
        print('   testing with ' + s +'.seq')
        seq = tom.load(s + '.seq')
        f = bz2.BZ2File(s + '.stree.bz2')
        streeForm = pickle.load(f)
        f.close()
        symbolSize = 2 if seq.nU() > 1 else 1
        # create SuffixTree in several random steps:
        sizes = sorted([rnd.randint(1, seq.length()-1) for i in range(5)]) + [seq.length()]
        stree = tom.STree(seq, symbolSize, symbolSize)
        for s in sizes:
            stree.extendTo(symbolSize * s)
            if not verifySuffixTreeCounts(seq.sub(0,symbolSize * s), stree):
                print('!!! SuffixTree test NOT passed !!! Substring count discrepancy.')
                ok = False
        streeFormNew = canonicalizeSuffixTree(stree)
        if not streeFormNew == streeForm:
            print('!!! SuffixTree test NOT passed !!! Something has changed.')
            ok = False
    if ok:
        print('...SuffixTree tests passed.')

#testSuffixTree()
