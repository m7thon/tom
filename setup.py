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
#
# obsolete options that were once required (older swig + older clang):
# cxx_flags += ['-stdlib=libc++']

try: from setuptools import setup,Extension,find_packages
except ImportError: from distutils.core import setup,Extension
import sys, os
import numpy as np

swig_flags = ['-c++', '-naturalvar', '-builtin', '-O', '-Iswig', '-DTOM_DEBUG']
if sys.version_info >= (3,): swig_flags += ['-py3']

def is_source_file(filename):
    extensions = ['.h', '.cpp', '.i', '.hpp']
    for ext in extensions:
        if filename.endswith(ext): return True
    return False

def make_depends():
    depends = []
    source_dirs = ['include/tom', 'swig']
    for source_dir in source_dirs:
        for root, dirs, files in os.walk(source_dir):
            for file in files:
                if is_source_file(file):
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
    ext_modules = [Extension('tom._tomlib', sources = ['swig/tomlib.i'], swig_opts = swig_flags,
                             language='c++', extra_compile_args = cxx_flags,
                             include_dirs = [np.get_include(), 'include/tom', 'include/external'],
                             depends = make_depends() )],
    packages = ['tom'],
    test_suite = 'tests',
)
