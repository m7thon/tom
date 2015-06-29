/* -*- C -*-  (not really, but good for syntax highlighting) */
%module(directors="1") tomlib
%naturalvar;
%feature("autodoc","1");

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
%include "std_vector.i"
%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;

%include "attribute.i"



%include "eigen3.i"
%init %{
    import_array();
%}

%include "tomdoc.i"

%include "../include/tom/tom.h"
%include "../include/tom/Macros.h"
%include "../include/tom/CerealTom.h"

%include "../include/tom/LinearAlgebra.h"
%include "../include/tom/Random.h"
//%shared_ptr(tom::Random)

%include "../include/tom/Sequence.h"
%shared_ptr(std::vector<tom::Sequence>);
//%shared_ptr(tom::Sequences);
%template(Sequences) std::vector<tom::Sequence>;

%include "../include/tom/stree/stree.h"
%shared_ptr(stree::STree);
%shared_ptr(std::vector<stree::Nidx>);
%template(NidxVector) std::vector<stree::Nidx>;

%include "../include/tom/stree/STreeCore.h"
%include "../include/tom/stree/STreeNode.h"
%include "../include/tom/stree/STreeIterators.h"

%include "../include/tom/PomdpTools.h"

%shared_ptr(tom::Oom);
%include "../include/tom/Oom.h"

%shared_ptr(tom::Hmm);
%include "../include/tom/Hmm.h"

%include "../include/tom/CoreSequences.h"
%include "../include/tom/Estimator.h"
%include "../include/tom/EfficiencySharpening.h"
