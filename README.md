Introduction
------------

Observable operator models are a new class of models for (controlled) stochastic processes that generalize HMMs / POMDPs and were developed by Herbert Jaeger.

This "Toolkit for observable operator modeling" (*tom*) aims to provide a reference implementation of the OOM methods developed in the [MINDS](minds.jacobs-university.de) research group of Jacobs University Bremen. Some benchmark problems as well as demo scripts will be included.

The core functionality is written in C++ for maximum performance. A Python interface to the library is included.

Requirements
------------

This toolkit relies on the following software:

- Eigen3: This is a convenient C++ matrix library that provides the basic
    linear algebra routines. Currently, the most recent (development) verision is required, which
    should correspond to a version >= 3.3.0.
- A recent C++ compiler supporting the current C++11 standard, e.g., gcc > 4.8, clang, ...
- SWIG: For generating the Python wrappers to the C++ code. A current version is required.
- Python and SciPy
- doxygen: To generate documentation from the source code

Furthermore, the toolkit makes use of the following open source software that is included for convenience:

- cereal (1.1.2): A C++ object serialization library. A modified version is included that uses a more current version
    of rapidjson and a different json writer object to produce more concise json output.
- rapidjson (1.0.2): A fast C++ json parser

Installation
------------

1. Install the required dependencies
2. Modify setup.cfg to set the path to the eigen3 header files if not in a standard path
3. run `python setup.py install`

Optionally,
- run `doxygen doc/tom.doxyfile` to generate the documentation in the subdirectory doc/html
- run `python setup.py clean --all` to clean up the build directory

Notes:
- The compiler and compile options can be modified by setting the environment variables `CC`, `CXX` and
  `CPPFLAGS`. Note that setup.py uses the compiler specified in `CC` for compiling and `CXX` for linking
  swig wrapped C++ extensions for unknown reasons, so always specify both. Example:  
  `CC=g++ CXX=g++ CPPFLAGS=-I/usr/local/include/eigen3 python setup.py install`.
- Useful flags / variables to define are:
  + `-DTOM_DEBUG` -- enable all debugging checks (asserts and bounds checking)
  + `-DTOM_NCHECK` -- disable bounds checking and such (enabled by default)
- You can also install tom into a user-local python package directory using  
  `python setup.py install --user`. For details, see the
  [distutils documentation](https://docs.python.org/3/install/index.html#alternate-installation).
- A Makefile is provided with some useful targets:
  + `make doc` -- make the documentation
  + `make clean` -- clean up
  + `make build` -- call `python setup.py build_ext`
  + `make install` -- call `python setup.py install`
  + `make debug_build_clang` -- compile with flags for debugging (`-DTOM_DEBUG`)
  + `make optimized_build_gcc5_on_OSX` -- compile the fastest version on my system with no bounds checking (`-DTOM_NCHECK`)
  + `make list` -- list all available targets
- Also, the Makefile allows building the toolkit without python's setuptools. For this
  1. Modify the first section in the Makefile to suit your local needs
  2. run `make build_without_setuptools`
  3. From python, include the directory python/tom into the search path:
     ```python
	import sys, os
	sys.path.append(os.path.abspath('.')+"/python/tom")
	import tom
    ```

Using the toolkit
-----------------

from python:
```python
import tom
```

Authors
-------

Tom is being developed by Michael Thon as part of his PhD thesis.

License
-------

Tom is licensed under the permissive open source [MIT License](http://opensource.org/licenses/MIT).
