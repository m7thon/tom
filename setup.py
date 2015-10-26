# ----------------------------------------------------------------------------
# Note: Instead of editing this file, you can also specify compiler options
#       from the command line as follows:
#
# You can specify the compiler to use by setting the environment variables
# `CC` and `CXX`, e.g., by invoking setup.py using
# CC=g++ CXX=g++ python setup.py install
#
# Similarly, you can set additional compile flags by setting the environment
# variables `CFLAGS` and `CXXFLAGS`, or simply `CPPFLAGS`, e.g., using
# CPPFLAGS="-march=native -O3" python setup.py install
# -----------------------------------------------------------------------------

cxx_flags = ['-std=c++11', '-Wno-unused-variable', '-Wno-sign-compare']

# optimization options (simply uncomment if this causes problems)
#cxx_flags += ['-march=native', '-O3']
#
# add minimal debug info on clang (gcc does not understand this)
#cxx_flags += ['-gline-tables-only']
#cxx_flags += ['-g0']
#
# when compiling on OSX with gcc and using -march=native (enabling avx instructions), use 
# the clang assembler, since the provided GNU as command does not support AVX instructions:
#cxx_flags += ['-Wa,-q']

try:
    from setuptools import setup,Extension
except ImportError:
    from distutils.core import setup,Extension

import os
import numpy as np

def make_depends():
    depends = []
    for source_dir in ['include', 'swig']:
        for root, dirs, files in os.walk(source_dir):
            for file in files:
                if os.path.splitext(file)[1] in ['.h', '.cpp', '.i', '.hpp']:
                    depends.append(os.path.join(root,file))
    return depends

setup(
    name = "tom",
    version = "0.4.0",
    author = "Michael Thon",
    author_email = "m7.thon@gmail.com",
    description = ("Toolkit for observable operator modeling"),
    license = "MIT",
    url = "https://gitlab.com/m7.thon/tom",
    package_dir = {'': 'python'},
    packages = ['tom', 'tom.hmm'],
    ext_modules = [Extension('tom._tomlib', sources = ['swig/tomlib_wrap.cpp'], language='c++',
                             extra_compile_args = cxx_flags,
                             include_dirs = [np.get_include(), 'include/tom', 'include/external'],
                             depends = make_depends() )],
    test_suite = 'tests',
)
