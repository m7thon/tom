import unittest
import tom

class PlainPySequence:
    """Simple list of symbols: [o_0, ..., o_{N-1}]."""

    def __init__(self, data):
        if type(data) is list:
            self.data = data
        else:
            self.data = list(range(data))

    def __repr__(self):
        return str(self.data)

    def length(self):
        return len(self.data)

    def raw(self):
        return self

    def equals(self, tom_seq):
        if tom_seq.rawSize() != len(self.data):
            return False
        for idx in range(len(self.data)):
            if self.data[idx] != tom_seq.rawAt(idx):
                return False
        return True

    def o(self, idx):
        return self.data[idx]

    def u(self, idx):
        return 0

    def sub(self, idx, length):
        return PlainPySequence(self.data[idx:idx + length:1 if length > 0 else -1])

    rawSub = sub

    def slice(self, begin, end, forwards):
        if begin == tom.NoIndex:
            begin = None
        if end == tom.NoIndex:
            end = None
        return PlainPySequence(self.data[begin:end:1 if forwards else -1])

    rawSlice = slice

    def reverse(self):
        return PlainPySequence(list(reversed(self.data)))

    def io(self, rev):
        if self.length() == 0:
            ret = IoPySequence([])
            ret.reversed = rev
            return ret
        dat = self.data[:]
        if (not rev and dat[0] > 0) or (rev and dat[0] < 0):
            dat = [None] + dat
        if len(dat) % 2 != 0:
            dat = dat + [None]
        i = dat[::2]
        o = dat[1::2]
        if rev:
            i, o = o, i
        ret = IoPySequence(list(zip(i, o)))
        ret.reversed = rev
        return ret


class IoPySequence:
    """ A list of io-pairs: [(u_0, o_0), ..., (u_{N-1}, o_{N-1})]. Note that for testing inputs will be negative
    and outputs positive!"""

    def __init__(self, data):
        if type(data) is list:
            self.data = data
        else:
            self.data = list(zip(range(-1, -(data + 1), -1), range(1, data + 1, 1)))
        self.reversed = False

    def __repr__(self):
        return str(self.data)

    def raw(self):
        if not self.reversed:
            return PlainPySequence([i for pair in self.data for i in pair if i is not None])
        else:
            return PlainPySequence([i for pair in self.data for i in reversed(pair) if i is not None])

    def length(self):
        return len(self.data)

    def equals(self, tom_seq):
        if not tom_seq.isIO():
            return False
        if tom_seq.length() != len(self.data):
            return False
        if len(self.data) > 0 and tom_seq.isReversed() != self.reversed:
            return False
        return self.raw().equals(tom_seq)

    def o(self, idx):
        return self.data[idx][1]

    def u(self, idx):
        return self.data[idx][0]

    def reverse(self):
        rev = IoPySequence(list(reversed(self.data)))
        rev.reversed = not self.reversed
        return rev

    def sub(self, idx, length):
        ret = IoPySequence(self.data[idx:idx + length:1 if length > 0 else -1])
        if length < 0:
            ret.reversed = not self.reversed
        else:
            ret.reversed = self.reversed
        return ret

    def slice(self, begin, end, forwards):
        if begin == tom.NoIndex:
            begin = None
        if end == tom.NoIndex:
            end = None
        ret = IoPySequence(self.data[begin:end:1 if forwards else -1])
        if not forwards:
            ret.reversed = not self.reversed
        else:
            ret.reversed = self.reversed
        return ret

    def rawSub(self, idx, length):
        return self.raw().sub(idx, length).io((self.reversed) == (length >= 0))

    def rawSlice(self, begin, end, forwards):
        return self.raw().slice(begin, end, forwards).io(self.reversed == forwards)


class TestSequence(unittest.TestCase):
    def create_subs(self):
        pass

    def cases(self):
        size = 5
        l = list(range(1, size + 1))
        py_base = PlainPySequence(l)
        tom_base = tom.Sequence(l, 10, 0)
        for b in range(size):
            for e in range(b, size + 1):
                try:
                    tom_seq = tom_base.rawSlice(b, e, True)
                    py_seq = py_base.slice(b, e, True)
                except:
                    self.assertTrue(False, "Error creating test cases by rawSlice by " + str(tom_base) + " [%d:%d:%d]" % (b, e, True))
                self.assertTrue(py_seq.equals(tom_seq), "Error creating test cases by rawSlice for " + str(tom_seq) + " and " + str(py_seq))
                yield (tom_seq, py_seq)
                py_seq = py_seq.reverse()
                tom_seq.reverse()
                self.assertTrue(py_seq.equals(tom_seq), "Error creating test cases by rawSlice for " + str(tom_seq) + " and " + str(py_seq))
                yield (tom_seq, py_seq)
        i = list(range(-1, -size - 1, -1))
        o = list(range(1, size + 1, 1))
        l = list(zip(i, o))
        x = []
        for p in l:
            x.extend(p)
        py_base = IoPySequence(l)
        tom_base = tom.Sequence(x, 10, 1)
        for b in range(2 * size):
            for e in range(b, 2 * size + 1):
                tom_seq = tom_base.rawSlice(b, e, True)
                py_seq = py_base.rawSlice(b, e, True)
                self.assertTrue(py_seq.equals(tom_seq), "Error creating test cases by rawSlice for " + str(tom_seq) + " and " + str(py_seq))
                yield (tom_seq, py_seq)
                py_seq = py_seq.reverse()
                tom_seq.reverse()
                self.assertTrue(py_seq.equals(tom_seq), "Error creating test cases by rawSlice for " + str(tom_seq) + " and " + str(py_seq))
                yield (tom_seq, py_seq)

    def test_json_io(self):
        for tom_seq, py_seq in self.cases():
            seq = tom.Sequence(tom_seq.toJSON())
            self.assertTrue(seq == tom_seq, "to and from json gives non-equal sequence for" + str(tom_seq) + " and " + str(seq))
            self.assertTrue(seq.nU() == tom_seq.nU(), "alphabet changed: " + str(tom_seq) + " and " + str(seq))
            self.assertTrue(seq.nO() == tom_seq.nO(), "alphabet changed: " + str(tom_seq) + " and " + str(seq))
        json = """{"Type":"Sequence","nU":1,"nO":10,"data":[-1,1,-2,2],"size":4}"""
        tom_seq = tom.Sequence(json)
        py_seq = IoPySequence([(-1, 1), (-2, 2)])
        self.assertTrue(tom_seq.nU() == 1 and tom_seq.nO() == 10 and py_seq.equals(tom_seq), "Error reading simple json-string")

    def test_copy(self):
        for tom_seq, py_seq in self.cases():
            seq = tom_seq.copy()
            self.assertTrue(seq == tom_seq, ".copy() not equal:" + str(tom_seq) + " and " + str(seq))
            self.assertTrue(seq.nU() == tom_seq.nU(), "alphabet changed: " + str(tom_seq) + " and " + str(seq))
            self.assertTrue(seq.nO() == tom_seq.nO(), "alphabet changed: " + str(tom_seq) + " and " + str(seq))

    def test_accessors(self):
        for tom_seq, py_seq in self.cases():
            for idx in range(-py_seq.length(), py_seq.length()):
                try:
                    self.assertTrue(tom_seq.o(idx) == py_seq.o(idx), ".o(%d) not correct: " % idx + str(tom_seq) + " and " + str(py_seq))
                except:
                    if py_seq.o(idx) is not None:
                        self.assertTrue(False, ".o(%d) should be %d: " % (idx, py_seq.o(idx)) + str(tom_seq) + " and " + str(py_seq))
                try:
                    self.assertTrue(tom_seq.u(idx) == py_seq.u(idx), ".u(%d) not correct: " % idx + str(tom_seq) + " and " + str(py_seq))
                except:
                    if py_seq.u(idx) is not None:
                        self.assertTrue(False, ".u(%d) should be %d: " % (idx, py_seq.u(idx)) + str(tom_seq) + " and " + str(py_seq))
            for idx in range(-py_seq.raw().length(), py_seq.raw().length()):
                self.assertTrue(tom_seq.rawAt(idx) == py_seq.raw().o(idx), "Error with rawAt: " + str(tom_seq))
                self.assertTrue(tom_seq.rawAt(idx) == tom_seq[idx], "Error with python []: " + str(tom_seq))
            self.assertTrue(list(tom_seq) == py_seq.raw().data, "Error with python iterator access: " + str(tom_seq))

    def test_rawSub(self):
        for tom_seq, py_seq in self.cases():
            for idx in list(range(py_seq.raw().length())):
                for l in list(range(py_seq.raw().length()-idx)):
                    self.assertTrue(py_seq.rawSub(idx, l).equals(tom_seq.rawSub(idx, l)), "Sub error: " + str(tom_seq) + " [%d:%d]" % (idx, l))
                for l in list(range(-1, -idx-1, -1)):
                    self.assertTrue(py_seq.rawSub(idx, l).equals(tom_seq.rawSub(idx, l)), "Sub error: " + str(tom_seq) + " [%d:%d]" % (idx, l))

    def test_sub(self):
        for tom_seq, py_seq in self.cases():
            for idx in list(range(py_seq.length())):
                for l in list(range(py_seq.length()-idx)):
                    self.assertTrue(py_seq.sub(idx, l).equals(tom_seq.sub(idx, l)), "Sub error: " + str(tom_seq) + " [%d:%d]" % (idx, l))
                for l in list(range(-1, -idx-1, -1)):
                    self.assertTrue(py_seq.sub(idx, l).equals(tom_seq.sub(idx, l)), "Sub error: " + str(tom_seq) + " [%d:%d]" % (idx, l))

    def test_rawSlice(self):
        for tom_seq, py_seq in self.cases():
            for b in [tom.NoIndex] + list(range(py_seq.raw().length())):
                if b == tom.NoIndex:
                    es = [tom.NoIndex] + list(range(py_seq.raw().length() + 1))
                else:
                    es = [tom.NoIndex] + list(range(b, py_seq.raw().length()+1))
                for e in es:
                    self.assertTrue(py_seq.rawSlice(b, e, True).equals(tom_seq.rawSlice(b, e, True)), "Slicing error: " + str(tom_seq) + " [%d:%d:%d]" % (b, e, True))
                    self.assertTrue(tom_seq.rawSlice(b, e, True) == tom_seq[None if b == tom.NoIndex else b:None if e == tom.NoIndex else e],
                                    "Python-slicing error: " + str(tom_seq) + " [%d:%d:%d]" % (b, e, True))
            for b in [tom.NoIndex] + list(range(py_seq.raw().length())):
                if b == tom.NoIndex:
                    es = [tom.NoIndex] + list(range(py_seq.raw().length()))
                else:
                    es = [tom.NoIndex] + list(range(0, b+1))
                for e in es:
                    self.assertTrue(py_seq.rawSlice(b, e, False).equals(tom_seq.rawSlice(b, e, False)), "Slicing error: " + str(tom_seq) + " [%d:%d:%d]" % (b, e, False))
                    self.assertTrue(tom_seq.rawSlice(b, e, False) == tom_seq[None if b == tom.NoIndex else b:None if e == tom.NoIndex else e:-1],
                                    "Python-slicing error: " + str(tom_seq) + " [%d:%d:%d]" % (b, e, False))

    def test_slice(self):
        for tom_seq, py_seq in self.cases():
            for b in [tom.NoIndex] + list(range(py_seq.length())):
                if b == tom.NoIndex:
                    es = [tom.NoIndex] + list(range(py_seq.length() + 1))
                else:
                    es = [tom.NoIndex] + list(range(b, py_seq.length()+1))
                for e in es:
                    self.assertTrue(py_seq.slice(b, e, True).equals(tom_seq.slice(b, e, True)), "Slicing error: " + str(tom_seq) + " [%d:%d:%d]" % (b, e, True))
            for b in [tom.NoIndex] + list(range(py_seq.length())):
                if b == tom.NoIndex:
                    es = [tom.NoIndex] + list(range(py_seq.length()))
                else:
                    es = [tom.NoIndex] + list(range(0, b+1))
                for e in es:
                    self.assertTrue(py_seq.slice(b, e, False).equals(tom_seq.slice(b, e, False)), "Slicing error: " + str(tom_seq) + " [%d:%d:%d]" % (b, e, False))


