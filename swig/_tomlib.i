/* -*- text -*-  (not really, but good for syntax highlighting) */
%module(directors="1") _tomlib
%naturalvar;
//%feature("autodoc","0");

%{ 
#define SWIG_FILE_WITH_INIT
#include "../include/tom/tom.h"
// Turn off an annoying warning:
#pragma GCC diagnostic ignored "-Warray-bounds"
%}

%include "stdint.i"

%include "std_shared_ptr.i"
%include "std_string.i"
%include "std_iostream.i"
%include "std_vector.i" // this should provide swig::stop_iteration
%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;

%include "attribute.i"

%include "eigen3.i"
%init %{
    import_array();
%}

%include "tomdoc.i"

%include "../include/tom/tom.h"
%include "../include/tom/CerealTom.h"

%include "../include/tom/StopCondition.h"
%include "../include/tom/LinearAlgebra.h"
%include "../include/tom/Random.h"
//%shared_ptr(tom::Random)

%include "../include/tom/Sequence.h"
// note to self: don't touch the following 3 lines, and don't ask...
%shared_ptr(std::vector<tom::Sequence>);
%template(Sequences) std::vector<tom::Sequence>;
%shared_ptr(tom::Sequences);

%include "../include/tom/stree/stree.h"
%include "../include/tom/stree/STreeCore.h"
%include "../include/tom/stree/STreeNode.h"
%include "../include/tom/stree/STreeIterators.h"

%include "../include/tom/Policy.h"

%shared_ptr(tom::Oom);
%include "../include/tom/Oom.h"

%shared_ptr(tom::Hmm);
%include "../include/tom/Hmm.h"

%include "../include/tom/CoreSequences.h"
%include "../include/tom/Estimator.h"
%include "../include/tom/EfficiencySharpening.h"
