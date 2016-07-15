import unittest
import tom
import os

current_dir = os.path.abspath(os.path.dirname(__file__)) + '/'

def wordsAreEqual(X, Y):
    if len(X) != len(Y): return False
    for x,y in zip(X,Y):
        if x != y: return False
    return True

class TestWords(unittest.TestCase):
    def test_wordsFromData(self):
        testSeqs = ['10_10000', '4_3_10000', '4x3_10000', 'pc_10000']
        for s in testSeqs:
            print('.', end='', flush=True)
            seq = tom.load(current_dir + s + '.seq.bz2')
            stree = tom.STree(seq)
            rstree = tom.STree(seq.reverse())
            wordSettings = [(0,0,2**k,0) for k in range(4, 10)]
            wordSettings += [(2,5,2**k,0) for k in range(4, 10)]
            wordSettings += [(2,5,2**k,0) for k in range(4, 10)]
            for wS in wordSettings:
                XY = tom.wordsFromData(stree, *wS)
                tom.sortWords(XY)
                XYr = tom.wordsFromData(rstree, *wS)
                tom.reverseWords(XYr)
                tom.sortWords(XYr)
                self.assertTrue(wordsAreEqual(XY, XYr),
                                "wordsFromData gives different results on reversed input for " + str(wS))
                Y = tom.wordsFromData(stree, *wS, prefixUnique=True)
                tom.sortWords(Y)
                Yr = tom.wordsFromData(rstree, *wS, suffixUnique=True)
                tom.reverseWords(Yr)
                tom.sortWords(Yr)
                self.assertTrue(wordsAreEqual(Y, Yr),
                                "wordsFromData gives different results on reversed input for characteristic words for " + str(wS))
                Xr = tom.wordsFromData(rstree, *wS, prefixUnique=True)
                tom.reverseWords(Xr)
                tom.sortWords(Xr)
                X = tom.wordsFromData(stree, *wS, suffixUnique=True)
                tom.sortWords(X)
                self.assertTrue(wordsAreEqual(X, Xr),
                                "wordsFromData gives different results on reversed input for indicative words for " + str(wS))
        print(' ', end='', flush=True)