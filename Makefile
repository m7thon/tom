include Makefile.inc

INCLUDE := $(EIGEN_INCLUDE) -I./include/tom -I./include/external

SRC := Makefile Makefile.inc
CXX_SRC := $(find include/tom -type f -name "*.h" -o -name "*.cpp" -o -name "*.hpp")
SWIG_SRC := $(find swig -type f -name "*.i")
PY_SRC := $(find python/tom -type f -name "*.py")

python/tom/_tomlib.so: swig/tomlib_wrap.cxx $(SRC) $(CXX_SRC)
	$(CXX) $(CXXFLAGS) $(OPT) $(PY_INCLUDE) $(INCLUDE) -shared -o python/tom/_tomlib.so swig/tomlib_wrap.cxx $(PY_LDFLAGS) $(LDFLAGS) 

swig/tomlib_wrap.cxx: $(SRC) $(CXX_SRC) $(SWIG_SRC)
	swig $(SWIG_FLAGS) -Iswig -outdir python/tom swig/tomlib.i

.PHONY: swig
swig:
	swig $(SWIG_FLAGS) -Iswig -outdir python/tom swig/tomlib.i

.PHONY: install
install: python/tom/_tomlib.so
	cd python; $(PYTHON) setup.py install

.PHONY: local
local: python/tom/_tomlib.so
	cd python; $(PYTHON) setup.py install --prefix=~/local

.PHONY: doc
doc:
	doxygen doc/tom.doxyfile

clean:
	rm -f python/tom/*.pyc
	rm -f python/tom/_tomlib.so python/tom/tomlib.py
	rm -f swig/tomlib_wrap.cxx swig/tomlib_wrap.h
	rm -rf python/build

.PHONY: stree
stree:
	swig $(SWIG_FLAGS) -Iswig -outdir lib swig/stree.i
	$(CXX) $(CXXFLAGS) $(OPT) $(PY_INCLUDE) $(INCLUDE) -shared -o lib/_stree.so swig/stree_wrap.cxx $(PY_LDFLAGS) $(LDFLAGS)

.PHONY: swig-docstrings
swig-docstrings:
	doc/doxy2swig.py -n -q doc/xml/index.xml swig/tom_doc.i

test: sandbox/test.cpp $(SRCS)
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDE) -o sandbox/test sandbox/test.cpp $(LDFLAGS)
