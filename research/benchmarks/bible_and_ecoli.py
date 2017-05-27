import tom
import urllib.request
from io import BytesIO
from zipfile import ZipFile
import bz2; import re

import random
random.seed(123456789)

# Download 'bible.txt' and 'E.coli' from the large canterbury corpus:
with urllib.request.urlopen('http://corpus.canterbury.ac.nz/resources/large.zip') as response:
    with ZipFile(BytesIO(response.read())) as zipfile:
        bible = zipfile.read('bible.txt').decode('UTF-8')
        ecoli = zipfile.read('E.coli').decode('UTF-8')

nTest = 5; lTest = 10**4

# Preprocess 'bible.tex' and generate training and test data for BIBLE27
bible = re.split(r'[\.\?!]+', bible)
random.shuffle(bible)
for idx in range(len(bible)):
    bible[idx] = re.sub(r'[ \.\?!,;:\-\(\)\'\n]+', r' ', bible[idx] + ' ').upper().lstrip()
symbolTableB27 = {' ':26}
for c in range(26):
    symbolTableB27[chr(c+65)]=c
print('BIBLE27: %d sentences' %(len(bible)))

tt = nTest * ['']
print('  length of test sequences: ', end='')
for test in range(nTest):
    while len(tt[test]) + len(bible[-1]) <= lTest:
        tt[test] += bible.pop()
    tom.save(tom.util.stringToSequence(tt[test], symbolTableB27)[0], 'BIBLE27_test%d' % test + '.seq.bz2')
    print('%d ' %(len(tt[test])), end='')
print()
tn = ''.join(bible)
tom.save(tom.util.stringToSequence(tn, symbolTableB27)[0], 'BIBLE27_train0.seq.bz2')
print('  length of train sequence:', len(tn))

# Preprocess 'E.coli' data and generate training and test data for ECOLI
symbolTable = {'g': 0, 'a': 1, 'c': 2, 't': 3}
snippletLength = len(ecoli) // nTest
tn = ''; tt = nTest * ['']
print('ECOLI')
print('  length of train sequences: ', end='')
for t in range(nTest):
    tn += ecoli[t * snippletLength : (t+1) * snippletLength - lTest]
    tt[t] = ecoli[(t+1) * snippletLength - lTest : (t+1) * snippletLength]
    tom.save(tom.util.stringToSequence(tt[t], symbolTable)[0], 'ECOLI_test%d' % t + '.seq.bz2')
    print('%d ' %(len(tt[t])), end='')
print()
tom.save(tom.util.stringToSequence(tn, symbolTable)[0], 'ECOLI_train0.seq.bz2')
print('  length of train sequence:', len(tn))
