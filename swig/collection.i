/* -*- C -*-  (not really, but good for syntax highlighting) */

%include "exception.i"

%{
#include <cassert>
#include <stdexcept>
%}

%define %RANGE_ERROR
{ try { $action }
	catch(std::out_of_range) { SWIG_exception(SWIG_IndexError, "Index out of bounds"); }
}
%enddef // %RANGE_ERROR

%define %LENGTH
%feature("python:slot", "sq_length", functype="lenfunc") __len__;
long __len__() const { return $self->size(); }
%enddef // LENGTH

%define %GETITEM(ITEM_TYPE)
%feature("python:slot", "mp_subscript", functype="binaryfunc") __getitem__;
%exception __getitem__ %RANGE_ERROR;
ITEM_TYPE __getitem__(long i) {
	if (i < 0) i = i + $self->size();
	if ((i < 0) or (i >= $self->size()))
		throw std::out_of_range("Index out of bounds");
	return (*($self))[i];
}
%enddef // %GETITEM(ITEM_TYPE)

%define %SETITEM(ITEM_TYPE)
%feature("python:slot", "mp_ass_subscript", functype="objobjargproc") __setitem__;
%exception __setitem__ %RANGE_ERROR;
void __setitem__(long i, const ITEM_TYPE& val) {
	if (i < 0) i = i + $self->size();
	if ((i < 0) or (i >= $self->size()))
		throw std::out_of_range("Index out of bounds");
	(*($self))[i] = val;
}
%enddef // %SETITEM(ITEM_TYPE)

%define %COLLECTION(ITEM_TYPE)
%LENGTH
%GETITEM(ITEM_TYPE)
%SETITEM(ITEM_TYPE)
%enddef // %COLLECTION(ITEM_TYPE)

%define %COLLECTION_RO(ITEM_TYPE)
%LENGTH
%GETITEM(ITEM_TYPE)
%enddef // %COLLECTION(ITEM_TYPE)
