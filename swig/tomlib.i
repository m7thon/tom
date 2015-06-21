/* -*- C -*-  (not really, but good for syntax highlighting) */
%module(directors="1") tomlib
%naturalvar;
%feature("autodoc","1");

// %include "tomdoc.i"

%{ 
#define SWIG_FILE_WITH_INIT
#include "../include/tom/tom.h"
// Turn off an annoying warning (HACK!):
#pragma GCC diagnostic ignored "-Warray-bounds"
%} 

%include "eigen3.i"

%include "../include/tom/tom.h"

%include "std_shared_ptr.i"

 //%include "cpointer.i"

%include "std_string.i"

%include "std_pair.i"
%template(IntIntPair) std::pair<int,int>;

%include "attribute.i"

%include "std_vector.i"
%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;

%include "stdint.i"
%include "std_iostream.i"

%init %{
  import_array();
%}

%include "../include/tom/Macros.h"
%include "../include/tom/CerealTom.h"
 //%pointer_class(tom::Symbol, SymbolPointer);

%include "../include/tom/LinearAlgebra.h"
%include "../include/tom/Random.h"

%include "../include/tom/Sequence.h"
%shared_ptr(std::vector<tom::Sequence>);
%shared_ptr(tom::Sequences);

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

 //%include "../include/tom/Estimators.h"
%include "../include/tom/Estimator.h"
//%include "../include/tom/Learner.h"

%include "../include/tom/EfficiencySharpening.h"
