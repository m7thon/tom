#ifndef TOM_H
#define TOM_H

//MARK: Macros

/** This is the lowest and default checking level and amounts to bounds checking on interface function, i.e., to catch user error (as opposed to internal debugging assertions). Errors will throw exceptions that should automatically be passed to python. Can be overridden by defining `TOM_NCHECK`. To enable all internal debug assertions (which will crash if they fail), define `TOM_DEBUG`. Note that these defines override (`>`) each other as follows:

    - `TOM_DEBUG > TOM_NCHECK > TOM_CHECK`
    - `TOM_DEBUG > NDEBUG`
*/
#define TOM_CHECK
#ifdef TOM_DEBUG
    #undef NDEBUG
    #undef TOM_NCHECK
#endif
#ifdef TOM_NCHECK
    #undef TOM_CHECK
#endif
#ifdef TOM_CHECK
    #define CHECK(...) __VA_ARGS__
#else
    #define CHECK(...)
#endif

/** Define a simple macro `SWIGCODE(...)` that facilitates writing SWIG code directly into the C++ source. This is hidden when not parsing by SWIG. */
#ifdef SWIG
    #define SWIGCODE(...) __VA_ARGS__
    #define NOSWIG(...)
#else
    #define SWIGCODE(...)
    #define NOSWIG(...) __VA_ARGS__
#endif

//#define VERBOSE

#ifdef VERBOSE
    #define LOOP_PROGRESS(name, i, end) {                                                        \
        if (i % ((unsigned long)end / 100) == 0) {                                               \
            fprintf(stderr, "\r%s...%3d%%", name, int(double(i) / end * 100));                   \
            fflush(stderr);                                                                      \
        }}
    #define LOOP_DONE(name) { fprintf(stderr, "\r%s...100%%\n", name); fflush(stderr); }
#else
    #define LOOP_PROGRESS(name, i, end)
    #define LOOP_DONE(name)
#endif

//MARK: C++ headers
#include <limits>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <deque>
#include <stack>
#include <queue>
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include <memory> // shared_ptr
#include <cstdint>
#include <cassert>
#include <stdexcept>
#include <random>

//MARK: external libraries
#include "Eigen/Core"
#include "Eigen/Dense"
namespace tom {
    using namespace Eigen;
}

#ifdef SWIG
namespace Eigen {
/** Return the number of threads used by the eigen3 backend (if compiled with openmp) */
int nbThreads();

/** Set the number of threads to be used by the eigen3 backend (if compiled with openmp) */
void setNbThreads(int n);
}
#endif

#include "cereal/cereal.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/string.hpp"
#include "cereal/archives/json.hpp"
#include "CerealTom.h"

/* The following C1 to C7 and PY1 to PY7 macros and are used to tweak the automatic documentation generation by doxygen for Python and C++. In code, the C macros just place their arguments verbatim, and the PY macros are ignored. */
#define C1(a) a
#define C2(a,b) a,b
#define C3(a,b,c) a,b,c
#define C4(a,b,c,d) a,b,c,d
#define C5(a,b,c,d,e) a,b,c,d,e
#define C6(a,b,c,d,e,f) a,b,c,d,e,f
#define C7(a,b,c,d,e,f,g) a,b,c,d,e,f,g
#define PY1(a)
#define PY2(a,b)
#define PY3(a,b,c)
#define PY4(a,b,c,d)
#define PY5(a,b,c,d,e)
#define PY6(a,b,c,d,e,f)
#define PY7(a,b,c,d,e,f,g)

//MARK: tom headers
namespace tom {
    typedef int Symbol;
}
#include "Random.h"
#include "Sequence.h"
#include "stree/stree.h"

#include "StopCondition.h"
#include "LinearAlgebra.h"
#include "Policy.h"
#include "Hmm.h"
#include "Oom.h"
#include "CoreSequences.h"
#include "Estimator.h"
#include "EfficiencySharpening.h"
#include "Implementations.h"

#endif // TOM_H
