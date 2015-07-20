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
    #define TOM_CHECK
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
#include "Eigen/Dense"
namespace tom {
    using namespace Eigen;
}
#include "cereal/cereal.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/string.hpp"
#include "cereal/archives/json.hpp"
#include "CerealTom.h"

//MARK: tom headers
namespace tom {
    typedef int Symbol;
}
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

//MARK: tom implementation
#include "Oom.cpp"
#include "LinearAlgebra.cpp"

#endif // TOM_H
