#ifndef TOM_H_
#define TOM_H_

/*! \mainpage Tom: Toolkit for observable operator modeling
 *
 * \section intro Introduction
 * Observable operator models are a new class of models for (controlled) stochastic processes that generalize HMMs / POMDPs and were developed by Herbert Jaeger.
 * 
 * This "Toolkit for observable operator modeling" (TOM) aims to provide a reference implementation of the OOM methods developed in the MINDS research group of Jacobs University Bremen. Some benchmark problems as well as demo scripts will be included.
 * 
 * The tom core functionality is written in C++ for maximum performance. A Python interface to the library is included.
 *
 * \section require Requirements
 * This toolkit relies on the following software:
 * - eigen3: This is a convenient C++ matrix library that provides the basic
     linear algebra routines. It is currently required to use the most recent
     (development) verision, which means a version strictly greater than 3.2.2
 * - A recent C++ compiler supporting the current C++11 standard, e.g., gcc > 4.8, clang, ...
 * - GNU make for building the toolkit
 * - SWIG: For generating the Python wrappers to the C++ code. Version 2.0.12 or later is required.
 * - Python/SciPy for scripting
 * - doxygen: To generate documentation from the source code
 *
 * \section install Installation
 * 
 * 1. Install the required dependencies
 * 2. Modify the Makefile.inc to set the correct compiler and paths to the dependencies
 * 3. run "make"
 * 4. run "make doc" to generate the documentation in the subdirectory doc/html
 *
 * \section usage Using the Toolkit
 *
 * from python, add the ./tom/lib directory to the system path:
 *
 * \>\>\> import sys, os
 *
 * \>\>\> sys.path.append(<path to tom\> + "/tom/lib")
 *
 * and import the toolkit:
 *
 * \>\>\> import tom
 *
 * \section authors Authors
 * Tom is being developed by Michael Thon as part of his PhD thesis.
 */

#include <limits>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <deque>
#include <vector>
#include <cmath>
#include <algorithm>

#include <memory> // shared_ptr
// #include "boost/shared_ptr.hpp"
// #include <tr1/memory>

#define SHARED_PTR std::shared_ptr

#include "../external/Eigen/Core"

#ifndef TOM_SYMBOL
#define TOM_SYMBOL int
#endif
#ifndef TOM_REAL
#define TOM_REAL double
#endif

#define STREE_STRING_TYPE tom::Sequence

namespace tom {
using namespace Eigen;
typedef TOM_SYMBOL Symbol;
typedef TOM_REAL Real;

typedef MatrixXi imat;
typedef MatrixXd mat;
typedef VectorXi ivec;
typedef VectorXd vec;
typedef RowVectorXi irvec;
typedef RowVectorXd rvec;

} // namespace tom

#include "CerealTom.h"
#include "Macros.h"
#include "Random.h"
#include "Sequence.h"

#include "stree/stree.h"

#include "LinearAlgebra.h"
#include "PomdpTools.h"
#include "Oom.h"
#include "Hmm.h"
#include "CoreSequences.h"
#include "Estimator.h"
#include "EfficiencySharpening.h"

// Include the implementation also:
#include "Oom.cpp"
#include "LinearAlgebra.cpp"

#endif // TOM_H_
