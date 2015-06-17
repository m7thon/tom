/* -*- C -*-  (not really, but good for syntax highlighting) */
%module stree
%naturalvar;
%feature("autodoc","1");
%{
#include "../include/tom/stree/stree.h"
%} 

%include "std_string.i"
%include "stdint.i"
%include "collection.i"

%nspace stree;

%include "../include/tom/stree/stree.h"


%shared_ptr(stree::STree);
%shared_ptr(std::vector<stree::Nidx>);
%template(NidxVector) std::vector<stree::Nidx>;

%include "../include/tom/stree/STreeCore.h"
%include "../include/tom/stree/STreeNode.h"
%include "../include/tom/stree/STreeIterators.h"
