#####################################################################
# WARNING! This makefile is just a wrapper around the setuptools
#          build system for python packages. The recommended way
#          to build and install is to:
#       1) Modify (if required) setup.cfg and setup.py
#       2) python setup.py build_ext
#       3) python setup.py install [--user]
#
# However, this Makefile provides some useful targets:
# * oldstyle:
#       this allowed to build the toolkit without the setuptools.
#       for this (and only for this!) modify the following part
#       of this makefile to reflect your configuration
# * clean:
#	clean up (this also removes the swig wrappers)
# * doc:
#       create the source code documentation -> doc/html/index.html
# * list:
#	list all available targets. Read this file for details
#####################################################################

# Path to eigen3 library if not in some standard location
EIGEN3_INCLUDE = -I/opt/local/include/eigen3

#####################################################################
# The following configuration options take effect *ONLY* for
# make build_without_setuptools
#####################################################################

# System libraries if not in a default location
LDFLAGS += -L/opt/local/lib
#LDFLAGS += -L/usr/local/lib

# The compiler to use. If not specified, the system default will be used. Note that a
# recent compiler supporting C++11 is required
#CXX = g++
#CXX = clang++
#CXX = g++-mp-4.9

# c++11 is required!
CXXFLAGS += -std=c++11

# If using clang or on older OSX (and using c++11), one may have to use libc++ ?!
#CXXFLAGS += -stdlib=libc++

# On OSX, when building with gcc, one needs to link with clang to support avx instructions
#CXXFLAGS += -Wa,-q

# Specify if you are using a non-standard python:
#PYTHON = python2.7
#PYTHON_CONFIG = python2.7-config
PYTHON = python
PYTHON_CONFIG = python-config

# Optimization flags:
#OPT = -O2 -march=msse2 -fPIC
OPT = -O3 -march=native -DNDEBUG
#OPT = -O -march=native -DNDEBUG -DTOM_NCHECK

#####################################################################
# End of basic configuration !
# You should not need to touch anything below here.
#####################################################################

# The following should work generically...
PYTHON_VERSION = $(shell ${PYTHON} -c 'import sys; print(sys.version_info[0])')
PY_INCLUDE = $(shell ${PYTHON_CONFIG} --includes)
NUMPY_INCLUDE = -I$(shell ${PYTHON} -c 'import numpy; print(numpy.get_include())')
PY_INCLUDE += ${NUMPY_INCLUDE} 
PY_LDFLAGS = $(shell ${PYTHON_CONFIG} --ldflags)

SWIG_FLAGS += -c++ -python -naturalvar -builtin -O
ifeq ($(PYTHON_VERSION),3)
  SWIG_FLAGS += -py3
endif

INCLUDE := $(EIGEN3_INCLUDE) -I./include/tom -I./include/external

CXX_SRC := $(shell find include/tom -type f -name '*.h' -o -name '*.cpp' -o -name '*.hpp')
SWIG_SRC := $(shell find swig -type f -name "*.i")
PY_SRC := $(shell find python/tom -type f -name "*.py")

.PHONY: oldstyle
oldstyle: python/tom/_tomlib.so
	@echo "-----------------------------------------------------------"
	@echo "You should now be able to use the toolkit. From python, do:"
	@echo ">>> import sys, os"
	@echo ">>> sys.path.append(os.path.abspath('.')+\"/python/tom\")"
	@echo ">>> import tom"
	@echo "-----------------------------------------------------------"

.PHONY: build
build:
	CPPFLAGS=$(EIGEN3_INCLUDE) python setup.py build_ext

.PHONY: build_clang
build_clang:
	CC=clang++ CXX=clang++ CPPFLAGS="$(EIGEN3_INCLUDE) -gline-tables-only -march=native -O3" python setup.py build_ext

.PHONY: build_gcc
build_gcc:
	CC=g++ CXX=g++ CPPFLAGS="$(EIGEN3_INCLUDE) -g0 -march=native -O3" python setup.py build_ext

.PHONY: build_gcc5_on_OSX
build_gcc5_on_OSX:
	CC=g++-mp-5 CXX=g++-mp-5 CPPFLAGS="$(EIGEN3_INCLUDE) -Wa,-q -Wa,-w -g0 -march=native -O3" python setup.py build_ext

.PHONY: deploy_fast
deploy_fast: swig/tomlib_wrap.cpp
	CC=clang++ CXX=clang++ CPPFLAGS="$(EIGEN3_INCLUDE) -gline-tables-only -O0 -DTOM_DEBUG" python setup.py install --user

.PHONY: deploy
deploy: swig/tomlib_wrap.cpp
	CC=clang++ CXX=clang++ CPPFLAGS="$(EIGEN3_INCLUDE) -gline-tables-only -march=native -O3" python setup.py install --user

.PHONY: test
test: swig/tomlib_wrap.cpp
	CC=clang++ CXX=clang++ CPPFLAGS="$(EIGEN3_INCLUDE) -gline-tables-only -march=native -O3" python setup.py install --user
	CC=clang++ CXX=clang++ CPPFLAGS="$(EIGEN3_INCLUDE) -gline-tables-only -march=native -O3" python setup.py test

.PHONY: install
install:
	python setup.py install

.PHONY: install_user
install_user:
	python setup.py install --user

.PHONY: optimized_build_clang
optimized_build_clang:
	CC=clang++ CXX=clang++ CPPFLAGS="$(EIGEN3_INCLUDE) -gline-tables-only -march=native -O3 -DTOM_NCHECK" python setup.py build_ext

.PHONY: optimized_build_gcc
optimized_build_gcc:
	CC=g++ CXX=g++ CPPFLAGS="$(EIGEN3_INCLUDE) -g0 -march=native -O3 -DTOM_NCHECK" python setup.py build_ext

.PHONY: optimized_build_gcc5_on_OSX
optimized_build_gcc5_on_OSX:
	CC=g++-mp-5 CXX=g++-mp-5 CPPFLAGS="$(EIGEN3_INCLUDE) -Wa,-q -Wa,-w -g0 -march=native -O3 -DTOM_NCHECK" python setup.py build_ext

.PHONY: docs
docs:
	doxygen doc/tom.doxyfile
	$(PYTHON) doc/doxy2swig/doxy2swig.py -focaq doc/xml/index.xml swig/tomdoc.i

.PHONY: doc_debug
doc_debug:
	doxygen doc/tom.doxdebug

.PHONY: swig
swig:
	swig $(SWIG_FLAGS) -Iswig -outdir python/tom -o swig/tomlib_wrap.cpp swig/tomlib.i

.PHONY: clean
clean:
	@rm -rf doc/html doc/xml
	@rm -f swig/tomdoc.i
	@rm -f swig/tomlib_wrap.cpp swig/tomlib_wrap.h swig/tomlib.py
	@rm -f python/tom/tomlib.py python/tom/_tomlib.so
	@$(PYTHON) setup.py -q clean --all 2>/dev/null
	@find . | grep -E "(__pycache__|\.pyc|\.egg-info)$$" | xargs rm -rf
	@rm -rf dist build
	@echo "Squeaky clean!"

doc/xml/index.xml: $(CXX_SRC)
	doxygen doc/tom.doxyfile

swig/tomdoc.i: doc/doxy2swig/doxy2swig.py doc/xml/index.xml
	$(PYTHON) doc/doxy2swig/doxy2swig.py -focaq doc/xml/index.xml swig/tomdoc.i

swig/tomlib_wrap.cpp: Makefile $(CXX_SRC) $(SWIG_SRC) swig/tomdoc.i
	swig $(SWIG_FLAGS) -Iswig -DTOM_DEBUG -outdir python/tom -o swig/tomlib_wrap.cpp swig/tomlib.i

python/tom/_tomlib.so: Makefile swig/tomlib_wrap.cpp
	$(CXX) $(CXXFLAGS) $(OPT) $(PY_INCLUDE) $(INCLUDE) -shared -o python/tom/_tomlib.so swig/tomlib_wrap.cpp $(PY_LDFLAGS) $(LDFLAGS)

.PHONY: sandbox
sandbox: sandbox/test.cpp $(SRCS)
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDE) -o sandbox/test sandbox/test.cpp $(LDFLAGS)

.PHONY: list
list:
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]" && $$3 ~ "^#  Phony") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$'
