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
- GNU make for building the toolkit
- SWIG: For generating the Python wrappers to the C++ code. A current version is required.
- Python/SciPy for scripting
- doxygen: To generate documentation from the source code

Furthermore, the toolkit makes use of the following open source software that is included for convenience:

- cereal (1.1.2): A C++ object serialization library. A modified version is included that uses a more current version
    of rapidjson and a different json writer object to produce more concise json output.
- rapidjson (1.0.2): A fast C++ json parser

Installation
------------

1. Install the required dependencies
2. Modify the Makefile.inc to set the correct compiler and paths to the dependencies  
   The default settings in Makefile.inc way work for you, though.
3. run `make`
4. run `make doc` to generate the documentation in the subdirectory doc/html

Using the toolkit
-----------------

from python, add the ./python/tom directory to the system path:
```python
import sys, os
sys.path.append(<path to tom> + '/python/tom')
import tom
```

Authors
-------

Tom is being developed by Michael Thon as part of his PhD thesis.

License
-------

Tom is licensed under the permissive open source [MIT License](http://opensource.org/licenses/MIT).
