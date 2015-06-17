#####################################################################
# Basic configuration for building tom. Please uncomment as required:
#####################################################################

# The compiler to use. If not specified, the system default will be used. Note that a
# recent compiler supporting C++11 is required
#CXX = g++
CXX = clang++
#CXX = g++-mp-4.9

# c++11 is now required!
CXXFLAGS += -std=c++11

# If using clang or on older OSX (and using c++11), one may have to use libc++ ?!
CXXFLAGS += -stdlib=libc++

# On OSX, when building with gcc, one needt to link with clang to support avx instructions
#CXXFLAGS += -Wa,-q

# Specify if you are using a non-standard python:
#PYTHON = python2.7
#PYTHON_CONFIG = python2.7-config
PYTHON = python
PYTHON_CONFIG = python-config

# The following should work generically...
PYTHON_VERSION = $(shell ${PYTHON} -c 'import sys; print(sys.version_info[0])')
PY_INCLUDE = $(shell ${PYTHON_CONFIG} --includes)
NUMPY_INCLUDE = -I$(shell ${PYTHON} -c 'import numpy; print(numpy.get_include())')
PY_INCLUDE += ${NUMPY_INCLUDE} 
PY_LDFLAGS = $(shell ${PYTHON_CONFIG} --ldflags)

# Path to eigen3 library if not in the external directory
EIGEN_INCLUDE = -I/opt/local/include/eigen3
EIGEN_INCLUDE += -I/usr/local/include/eigen3

# System libraries if not in a default location
LDFLAGS += -L/opt/local/lib
LDFLAGS += -L/usr/local/lib

# Optimization flags:
OPT = -O3 -march=native -DNDEBUG # -fPIC
#OPT = -O2 -march=native -fPIC

#####################################################################
# End of basic configuration !
#####################################################################

SWIG_FLAGS += -c++ -python -naturalvar -builtin -O -modern -modernargs -directors
ifeq ($(PYTHON_VERSION),3)
  SWIG_FLAGS += -py3
endif
