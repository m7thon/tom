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
    packages = ['tom', 'tom.hmm', 'tom.learn', 'tom.linalg', 'tom.sequence', 'tom.util', 'tom.io'],
    ext_modules = [Extension('tom.__tomlib', sources = ['swig/_tomlib_wrap.cpp'], language='c++',
                             extra_compile_args = ['-std=c++11', '-Wno-unused-variable', '-Wno-sign-compare'],
                             include_dirs = [np.get_include(), 'include/tom', 'include/external'],
                             depends = make_depends() )],
    test_suite = 'tests',
)
