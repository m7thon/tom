#####################################################################
# WARNING! This makefile is just a wrapper around the setuptools
#          build system for python packages. The recommended way
#          to build and install is to:
#       1) Modify (if required) setup.cfg and setup.py
#       2) python setup.py build_ext
#       3) python setup.py install [--user]
#
# However, this Makefile provides some useful targets:
# * build_without_setuptools:
#       this allowd to build the toolkit without the setuptools.
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
EIGEN3_INCLUDE = -I/usr/local/include/eigen3

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
OPT = -O3 -march=native -DNDEBUG -fPIC
#OPT = -O2 -march=native -fPIC

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

CXX_SRC := $(find include/tom -type f -name "*.h" -o -name "*.cpp" -o -name "*.hpp")
SWIG_SRC := $(find swig -type f -name "*.i")
PY_SRC := $(find python/tom -type f -name "*.py")

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
	CC=g++ CXX=g++ CPPFLAGS="$(EIGEN3_INCLUDE) -Wa,-q -g0 -march=native -O3" python setup.py build_ext

.PHONY: build_without_setuptools
build_without_setuptools: python/tom/_tomlib.so
	@echo "-----------------------------------------------------------"
	@echo "You should now be able to use the toolkit. From python, do:"
	@echo ">>> import sys, os"
	@echo ">>> sys.path.append(os.path.abspath('.')+\"/python/tom\")"
	@echo ">>> import tom"

.PHONY: debug_build_clang
debug_build_clang:
	CC=clang++ CXX=clang++ CPPFLAGS="$(EIGEN3_INCLUDE) -gline-tables-only -O0 -DTOM_DEBUG" python setup.py build_ext

.PHONY: install
install:
	python setup.py install

.PHONY: install_user
install_user:
	python setup.py install --user

.PHONY: swig
swig:
	swig $(SWIG_FLAGS) -Iswig -outdir python/tom swig/tomlib.i

.PHONY: doc
doc:
	doxygen doc/tom.doxyfile

.PHONY: clean
clean:
	rm -rf dist
	rm -f python/tom/*.pyc
	rm -f python/tom/_tomlib.so python/tom/tomlib.py
	rm -f swig/tomlib_wrap.cxx swig/tomlib_wrap.cpp swig/tomlib_wrap.h
	rm -rf build
	python setup.py clean --all

.PHONY: swig-docstrings
swig-docstrings:
	doc/doxy2swig.py -n -q doc/xml/index.xml swig/tom_doc.i

# Not a target:
python/tom/_tomlib.so: Makefile swig/tomlib_wrap.cxx $(CXX_SRC)
	$(CXX) $(CXXFLAGS) $(OPT) $(PY_INCLUDE) $(INCLUDE) -shared -o python/tom/_tomlib.so swig/tomlib_wrap.cxx $(PY_LDFLAGS) $(LDFLAGS)

# Not a target:
swig/tomlib_wrap.cxx: $(CXX_SRC) $(SWIG_SRC)
	swig $(SWIG_FLAGS) -Iswig -DTOM_DEBUG -outdir python/tom swig/tomlib.i

# Not a target:
.PHONY: test
test: sandbox/test.cpp $(SRCS)
	$(CXX) $(CXXFLAGS) $(OPT) $(INCLUDE) -o sandbox/test sandbox/test.cpp $(LDFLAGS)

.PHONY: list
list:
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$'
