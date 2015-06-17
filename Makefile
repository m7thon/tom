include Makefile.inc

SRCs = Makefile Makefile.inc
SRCS += include/tom/stree/stree.h include/tom/stree/RBTree.h include/tom/stree/STreeCore.h include/tom/stree/STreeCore.cpp include/tom/stree/STreeNode.h include/tom/stree/STreeIterators.h
SRCS += include/tom/tom.h include/tom/Macros.h include/tom/CerealTom.h include/tom/Random.h include/tom/PomdpTools.h include/tom/LinearAlgebra.h include/tom/LinearAlgebra.cpp include/tom/Sequence.h include/tom/Oom.h include/tom/Hmm.h include/tom/Oom.cpp include/tom/Estimator.h include/tom/CoreSequences.h include/tom/EfficiencySharpening.h

PY_SRCS = python/tom.py python/tomio.py python/tomseqs.py python/tomlearn.py

SWIG_SRCS = swig/eigen3.i swig/collection.i swig/stree.i swig/tomlib.i

tom: $(PY_SRCS) lib/_tomlib.so
	cp python/*.py lib/.

.PHONY: swig
swig:
	swig $(SWIG_FLAGS) -Iswig -outdir lib swig/tomlib.i

swig/tomlib_wrap.cxx: $(SRCS) $(SWIG_SRCS)
	swig $(SWIG_FLAGS) -Iswig -outdir lib swig/tomlib.i

lib/_tomlib.so: $(SRCS) swig/tomlib_wrap.cxx $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPT) $(PY_INCLUDE) $(INCLUDE) -shared -o lib/_tomlib.so swig/tomlib_wrap.cxx $(PY_LDFLAGS) $(LDFLAGS)

.PHONY: doc
doc:
	doxygen doc/tom.doxyfile

clean:
	rm -f include/tom/*~ include/tom/stree/*~ swig/*~ python/*~
	rm -f swig/*.o
	rm -f lib/tomlib.py lib/tom.py lib/tomio.py lib/tomseqs.py lib/tomlearn.py lib/*.pyc
	rm -f lib/_tomlib.so lib/stree.py lib/_stree.so
	rm -f swig/tomlib_wrap.cxx swig/stree_wrap.cxx 

.PHONY: stree
stree:
	swig $(SWIG_FLAGS) -Iswig -outdir lib swig/stree.i
	$(CXX) $(CXXFLAGS) $(OPT) $(PY_INCLUDE) $(INCLUDE) -shared -o lib/_stree.so swig/stree_wrap.cxx $(PY_LDFLAGS) $(LDFLAGS)

.PHONY: swig-docstrings
swig-docstrings:
	doc/doxy2swig.py -n -q doc/xml/index.xml swig/tom_doc.i

test: sandbox/test.cpp $(SRCS)
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDE) -o sandbox/test sandbox/test.cpp $(LDFLAGS)
