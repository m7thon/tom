#ifndef TOM_H
#define TOM_H

#include "Macros.h"

#include <limits>
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

#include "Eigen/Core"


namespace tom {
using namespace Eigen;
typedef int Symbol;
typedef double Real;
} // namespace tom

#include "cereal/cereal.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/string.hpp"
#include "cereal/archives/json.hpp"
#include "CerealTom.h"

#include <random>
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
#include "Eigen/LU"
#include "Oom.cpp"

#include "Eigen/Cholesky"
#include "Eigen/SVD"
#include "Eigen/QR"
#include "Eigen/Eigenvalues"
#include "LinearAlgebra.cpp"

#endif // TOM_H
