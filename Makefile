PYTHON = python

SWIG_FLAGS = -c++ -python -naturalvar -builtin -O

CXX_SRC := $(shell find include/tom -type f -name '*.h' -o -name '*.cpp' -o -name '*.hpp')
SWIG_SRC := $(shell find swig -type f -name "*.i")
PY_SRC := $(shell find python/tom -type f -name "*.py")

.PHONY: build
.PHONY: install install_user
.PHONY: test debug
.PHONY: doc swig clean cleanall
.PHONY: deploy deploy_quickly deploy_optimized

build: swig/tomlib_wrap.cpp
	$(PYTHON) setup.py build_ext

install: swig/tomlib_wrap.cpp
	$(PYTHON) setup.py install

install_user: swig/tomlib_wrap.cpp
	$(PYTHON) setup.py install --user

test: swig/tomlib_wrap.cpp
	CPPFLAGS="-DTOM_DEBUG -g" $(PYTHON) setup.py test

debug: swig/tomlib_wrap.cpp
	CPPFLAGS="-DTOM_DEBUG -g" $(PYTHON) setup.py build_ext

deploy: swig/tomlib_wrap.cpp
	CC=clang++ CXX=clang++ CPPFLAGS="-g -march=native -O2 -DTOM_DEBUG" $(PYTHON) setup.py install --user

deploy_quickly: swig/tomlib_wrap.cpp
	CC=clang++ CXX=clang++ CPPFLAGS="-gline-tables-only -O0 -DTOM_DEBUG" $(PYTHON) setup.py install --user

deploy_optimized: swig/tomlib_wrap.cpp
	CC=g++-mp-5 CXX=g++-mp-5 CPPFLAGS="-Wa,-q -Wa,-w -g0 -march=native -O3 -DTOM_NCHECK" $(PYTHON) setup.py install --user

doc:
	doxygen doc/tom.doxyfile
	$(PYTHON) doc/doxy2swig/doxy2swig.py -focaq doc/xml/index.xml swig/tomdoc.i

swig: swig/tomdoc.i
	swig $(SWIG_FLAGS) -Iswig -outdir python/tom -o swig/tomlib_wrap.cpp swig/tomlib.i

clean:
	@rm -f python/tom/_tomlib.so
	@$(PYTHON) setup.py -q clean --all 2>/dev/null
	@find . | grep -E "(__pycache__|\.pyc|\.egg-info)$$" | xargs rm -rf
	@rm -rf dist build

cleanall: clean
	@rm -rf doc/html doc/xml
	@rm -f swig/tomdoc.i
	@rm -f swig/tomlib_wrap.cpp swig/tomlib_wrap.h python/tom/tomlib.py

doc/xml/index.xml: $(CXX_SRC)
	doxygen doc/tom.doxyfile

swig/tomdoc.i: doc/doxy2swig/doxy2swig.py doc/xml/index.xml
	$(PYTHON) doc/doxy2swig/doxy2swig.py -focaq doc/xml/index.xml swig/tomdoc.i

swig/tomlib_wrap.cpp: Makefile $(CXX_SRC) $(SWIG_SRC) swig/tomdoc.i
	swig $(SWIG_FLAGS) -Iswig -outdir python/tom -o swig/tomlib_wrap.cpp swig/tomlib.i
