/* -*- C++ -*-  (not really, but good for syntax highlighting) */
#ifdef SWIGPYTHON
%naturalvar;

%{
#include "Eigen/Core"
%}

// numpy.i is from https://github.com/numpy/numpy/tree/master/doc/swig
%include "numpy.i"


/*****************************************************************************
  Utility functions used for Eigen - NumPy conversions:
 ****************************************************************************/
%fragment("Eigen_NumPy_Utilities", "header", fragment="NumPy_Fragments")
{
/* Functions to handle Eigen objects wrapped inside a PyCapsule or PyCObject */
%#if NPY_API_VERSION < 0x00000007
  %#define array_setbase(a,b)     (PyArray_BASE(a)=b)
%#else
  %#define array_setbase(a,b)     (PyArray_SetBaseObject((PyArrayObject *)a,b))
%#endif

%#ifdef SWIGPY_USE_CAPSULE
  %#define encapsulate(cobj, ...) (PyCapsule_New(cobj,SWIGPY_CAPSULE_NAME, __VA_ARGS__))
	template <typename T_Ptr>
	void clean(PyObject* obj) {
		void* data = (void*) PyCapsule_GetPointer(obj,SWIGPY_CAPSULE_NAME);
		if (data != NULL) delete (T_Ptr) data;
	}
%#else
  %#define encapsulate(cobj, ...) (PyCObject_FromVoidPtr(cobj, __VA_ARGS__))
	template <typename T_Ptr>
	void clean(void* obj) { if (data != NULL) delete (T_Ptr) data; }
%#endif
}

/*****************************************************************************
  Some convenience typemaps:
 ****************************************************************************/
%{
	typedef Eigen::Map<Eigen::Matrix<    int, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  MatrixMi;
	typedef Eigen::Map<Eigen::Matrix<   long, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  MatrixMl;
	typedef Eigen::Map<Eigen::Matrix<  float, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  MatrixMf;
	typedef Eigen::Map<Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  MatrixMd;
	typedef Eigen::Map<Eigen::Array<     int, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  ArrayMMi;
	typedef Eigen::Map<Eigen::Array<    long, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  ArrayMMl;
	typedef Eigen::Map<Eigen::Array<   float, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  ArrayMMf;
	typedef Eigen::Map<Eigen::Array<  double, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  ArrayMMd;
	typedef Eigen::Map<Eigen::Array<     int, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> >  ArrayMi;
	typedef Eigen::Map<Eigen::Array<    long, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> >  ArrayMl;
	typedef Eigen::Map<Eigen::Array<   float, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> >  ArrayMf;
	typedef Eigen::Map<Eigen::Array<  double, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> >  ArrayMd;
%}

typedef Eigen::Map<Eigen::Matrix<    int, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  MatrixMi;
typedef Eigen::Map<Eigen::Matrix<   long, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  MatrixMl;
typedef Eigen::Map<Eigen::Matrix<  float, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  MatrixMf;
typedef Eigen::Map<Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  MatrixMd;
typedef Eigen::Map<Eigen::Array<     int, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  ArrayMMi;
typedef Eigen::Map<Eigen::Array<    long, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  ArrayMMl;
typedef Eigen::Map<Eigen::Array<   float, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  ArrayMMf;
typedef Eigen::Map<Eigen::Array<  double, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >  ArrayMMd;
typedef Eigen::Map<Eigen::Array<     int, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> >  ArrayMi;
typedef Eigen::Map<Eigen::Array<    long, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> >  ArrayMl;
typedef Eigen::Map<Eigen::Array<   float, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> >  ArrayMf;
typedef Eigen::Map<Eigen::Array<  double, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> >  ArrayMd;

namespace Eigen {
%ignore Dynamic;
  class Dynamic;
  template<typename T, typename D1, typename D2> class Matrix;
  template<typename T, typename D1, typename D2> class Array;
  typedef Matrix<   int, Dynamic, Dynamic> MatrixXi;
  typedef Matrix< float, Dynamic, Dynamic> MatrixXf;
  typedef Matrix<double, Dynamic, Dynamic> MatrixXd;
  typedef Matrix<   int, Dynamic, 1> VectorXi;
  typedef Matrix< float, Dynamic, 1> VectorXf;
  typedef Matrix<double, Dynamic, 1> VectorXd;
  typedef Matrix<   int, 1, Dynamic> RowVectorXi;
  typedef Matrix< float, 1, Dynamic> RowVectorXf;
  typedef Matrix<double, 1, Dynamic> RowVectorXd;
  typedef Array<   int, Dynamic, Dynamic> ArrayXXi;
  typedef Array< float, Dynamic, Dynamic> ArrayXXf;
  typedef Array<double, Dynamic, Dynamic> ArrayXXd;
  typedef Array<   int, Dynamic, 1> ArrayXi;
  typedef Array< float, Dynamic, 1> ArrayXf;
  typedef Array<double, Dynamic, 1> ArrayXd;
}


/*****************************************************************************
  The actual typemaps follow:
 ****************************************************************************/
%define %eigen_numpy_typemaps(DATA_TYPE, DATA_TYPECODE, DATA_TYPEPRECEDENCE)

%typemap(typecheck, precedence = DATA_TYPEPRECEDENCE, fragment="Eigen_NumPy_Utilities")
const Eigen::MatrixBase<Eigen::Map<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::ArrayBase <Eigen::Map<Eigen::Array <DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::ArrayBase <Eigen::Map<Eigen::Array <DATA_TYPE, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> > >&,
const Eigen::DenseBase <Eigen::Map<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::DenseBase <Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::DenseBase <Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> > >&,
const Eigen::EigenBase <Eigen::Map<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::EigenBase <Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::EigenBase <Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> > >&,
const Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>&,
const Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>&
{
  $1 = is_array($input) && type_match(array_type($input), DATA_TYPECODE);
}


/******************************
  Input const & typemaps     // Passing a wrapper object (Eigen::Map), but without copying the data
/*****************************/
%typemap(in, fragment="Eigen_NumPy_Utilities")
const Eigen::MatrixBase<Eigen::Map<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::DenseBase <Eigen::Map<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::EigenBase <Eigen::Map<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&
{
	PyArrayObject * ary = obj_to_array_no_conversion($input, DATA_TYPECODE);
	if (ary == NULL) SWIG_fail;
	if (array_numdims(ary) != 2) { PyErr_SetString(PyExc_ValueError, "array must be 2-dimensional"); SWIG_fail; }
	int rows = array_size(ary, 0), cols = array_size(ary, 1);
	int outer = array_stride(ary,1)/PyArray_ITEMSIZE(ary), inner = array_stride(ary,0)/PyArray_ITEMSIZE(ary);
  $1 = new Eigen::Map<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > ((DATA_TYPE*) array_data(ary), rows, cols, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(outer, inner));
}
%typemap(in, fragment="Eigen_NumPy_Utilities")
const Eigen::ArrayBase <Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::DenseBase <Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::EigenBase <Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&
{
	PyArrayObject * ary = obj_to_array_no_conversion($input, DATA_TYPECODE);
	if (ary == NULL) SWIG_fail;
	if (array_numdims(ary) != 2) { PyErr_SetString(PyExc_ValueError, "array must be 2-dimensional"); SWIG_fail; }
	int rows = array_size(ary, 0), cols = array_size(ary, 1);
	int outer = array_stride(ary,1)/PyArray_ITEMSIZE(ary), inner = array_stride(ary,0)/PyArray_ITEMSIZE(ary);
  $1 = new Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > ((DATA_TYPE*) array_data(ary), rows, cols, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(outer, inner));
}
%typemap(in, fragment="Eigen_NumPy_Utilities")
const Eigen::ArrayBase <Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> > >&,
const Eigen::DenseBase <Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> > >&,
const Eigen::EigenBase <Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> > >&
{
	PyArrayObject * ary = obj_to_array_no_conversion($input, DATA_TYPECODE);
	if (ary == NULL) SWIG_fail;
	if (array_numdims(ary) != 1) { PyErr_SetString(PyExc_ValueError, "array must be 1-dimensional"); SWIG_fail; }
	int size = array_size(ary, 0);
	int inner = array_stride(ary,0)/PyArray_ITEMSIZE(ary);
  $1 = new Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> > ((DATA_TYPE*) array_data(ary), size, Eigen::InnerStride<Eigen::Dynamic>(inner));
}
%typemap(freearg)
const Eigen::MatrixBase<Eigen::Map<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::ArrayBase <Eigen::Map<Eigen::Array <DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::ArrayBase <Eigen::Map<Eigen::Array <DATA_TYPE, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> > >&,
const Eigen::DenseBase <Eigen::Map<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::DenseBase <Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::DenseBase <Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> > >&,
const Eigen::EigenBase <Eigen::Map<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::EigenBase <Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > >&,
const Eigen::EigenBase <Eigen::Map<Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>, 0, Eigen::InnerStride<Eigen::Dynamic> > >&
{ if ($1) delete $1; }


/******************************
  Input const & typemaps     // This involves a temporary, i.e., a full copy of the data, but passes a true Eigen::MatrixXt object instead of an Eigen::Map.
/*****************************/
%typemap(in, fragment="Eigen_NumPy_Utilities")
  const Eigen::MatrixBase<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic> >& (Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic> temp),
  const Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>& (Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic> temp),
  const Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, 1>& (Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, 1> temp),
  const Eigen::Matrix<DATA_TYPE, 1, Eigen::Dynamic>& (Eigen::Matrix<DATA_TYPE, 1, Eigen::Dynamic> temp),
  const Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>& (Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic> temp),
  const Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>& (Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1> temp)
{
	PyArrayObject * ary = obj_to_array_no_conversion($input, DATA_TYPECODE);
	if (ary == NULL) SWIG_fail;
	if(array_numdims(ary) != 2) { PyErr_SetString(PyExc_ValueError, "array must be 2-dimensional"); SWIG_fail; }
	int rows = array_size(ary, 0), cols = array_size(ary, 1);
	int outer = array_stride(ary,1)/PyArray_ITEMSIZE(ary), inner = array_stride(ary,0)/PyArray_ITEMSIZE(ary);
  temp = Eigen::Map<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > ((DATA_TYPE*) array_data(ary), rows, cols, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(outer, inner));
	$1 = &temp;
}

/******************************
  Input value typemaps       // This involves at least one temporary
/*****************************/
%typemap(in, optimal=1, fragment="Eigen_NumPy_Utilities")
  Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>,
  Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, 1>,
  Eigen::Matrix<DATA_TYPE, 1, Eigen::Dynamic>,
  Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>,
  Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>
{
	PyArrayObject * ary = obj_to_array_no_conversion($input, DATA_TYPECODE);
	if (ary == NULL) SWIG_fail;
	if(array_numdims(ary) != 2) { PyErr_SetString(PyExc_ValueError, "array must be 2-dimensional"); SWIG_fail; }
	int rows = array_size(ary, 0), cols = array_size(ary, 1);
	int outer = array_stride(ary,1)/PyArray_ITEMSIZE(ary), inner = array_stride(ary,0)/PyArray_ITEMSIZE(ary);
  $1 = Eigen::Map<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > ((DATA_TYPE*) array_data(ary), rows, cols, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(outer, inner));
}

/******************************
  Output & typemaps          // This does not pass ownership to Python
/*****************************/
%typemap(out, fragment="Eigen_NumPy_Utilities")
  Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>&,
  Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, 1>&,
	Eigen::Matrix<DATA_TYPE, 1, Eigen::Dynamic>&,
  Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>&,
  Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>&
{
	npy_intp dims[2] = { $1->rows(), $1->cols() };
	$result = PyArray_New(&PyArray_Type, 2, dims, DATA_TYPECODE, NULL, (void*) const_cast<$1_ltype>($1)->data(), 0, NPY_ARRAY_FARRAY, NULL);
	if (!$result) SWIG_fail;
  array_setbase($result,encapsulate($1, NULL));
}

/******************************
  Output (const) & typemaps  //
/*****************************/
%typemap(out, fragment="Eigen_NumPy_Utilities")
  const Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>&,
  const Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, 1>&,
	const Eigen::Matrix<DATA_TYPE, 1, Eigen::Dynamic>&,
  const Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>&,
  const Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>&
{
	npy_intp dims[2] = { $1->rows(), $1->cols() };
	$result = PyArray_New(&PyArray_Type, 2, dims, DATA_TYPECODE, NULL, (void*) const_cast<$1_ltype>($1)->data(), 0, NPY_ARRAY_FARRAY_RO, NULL);
	if (!$result) SWIG_fail;
  array_setbase($result,encapsulate($1, NULL));
}

/******************************
  Output * typemaps          // This passes ownership to Python
/*****************************/
%typemap(out, fragment="Eigen_NumPy_Utilities")
  Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>*,
  Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, 1>*,
	Eigen::Matrix<DATA_TYPE, 1, Eigen::Dynamic>*,
  Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>*,
  Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>*
{
	npy_intp dims[2] = { $1->rows(), $1->cols() };
	$result = PyArray_New(&PyArray_Type, 2, dims, DATA_TYPECODE, NULL, (void*) $1->data(), 0, NPY_ARRAY_FARRAY, NULL);
	if (!$result) SWIG_fail;
  array_setbase($result, encapsulate($1, clean<$1_ltype>));
}

/******************************
  Output SHARED_PTR typemaps //
/*****************************/
%typemap(out, optimal = "1", fragment="Eigen_NumPy_Utilities")
  SHARED_PTR<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic> >,
  SHARED_PTR<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, 1> >,
  SHARED_PTR<Eigen::Matrix<DATA_TYPE, 1, Eigen::Dynamic> >,
  SHARED_PTR<Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic> >,
  SHARED_PTR<Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1> >
{
  $1_ltype* newSharedPtr = new $1_ltype((const $1_ltype &)$1);
	npy_intp dims[2] = { $1->rows(), $1->cols() };
	$result = PyArray_New(&PyArray_Type, 2, dims, DATA_TYPECODE, NULL, (void*) $1->data(), 0, NPY_ARRAY_FARRAY, NULL);
	if (!$result) SWIG_fail;
  array_setbase($result, encapsulate(newSharedPtr, clean<$1_ltype*>));
}

/******************************
  Output value typemaps      // This creates a copy and passes ownership of the copy to Python
/*****************************/
%typemap(out, optimal = "1", fragment="Eigen_NumPy_Utilities")
  Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>,
  Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, 1>,
	Eigen::Matrix<DATA_TYPE, 1, Eigen::Dynamic>,
  Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>,
  Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>
{
  $1_ltype* temp = new $1_ltype((const $1_ltype &)$1);
	npy_intp dims[2] = { temp->rows(), temp->cols() };
	$result = PyArray_New(&PyArray_Type, 2, dims, DATA_TYPECODE, NULL, (void*) temp->data(), 0, NPY_ARRAY_FARRAY, NULL);
	if (!$result) { delete temp; SWIG_fail; }
  array_setbase($result, encapsulate(temp, clean<$1_ltype*>));
}

/******************************
  Argout & typemaps          //
/*****************************/
%typemap(in, numinputs=0)
  const Eigen::MatrixBase<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic> >& OUTPUT
{
	$1 = new Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>();
}
%typemap(argout)
  const Eigen::MatrixBase<Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic> >& OUTPUT
{
	npy_intp dims[2] = { $1->rows(), $1->cols() };
	PyObject* res = PyArray_New(&PyArray_Type, 2, dims, DATA_TYPECODE, NULL, (void*) $1->derived().data(), 0, NPY_ARRAY_FARRAY, NULL);
	if (!res) { delete $1; SWIG_fail; }
	array_setbase(res, encapsulate($1, clean<Eigen::Matrix<DATA_TYPE,Eigen::Dynamic,Eigen::Dynamic>*>));
	$result = SWIG_Python_AppendOutput($result,res);
}
%typemap(in, numinputs=0)
  Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>& OUTPUT,
  Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, 1>& OUTPUT,
	Eigen::Matrix<DATA_TYPE, 1, Eigen::Dynamic>& OUTPUT,
  Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>& OUTPUT,
  Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>& OUTPUT
{
	$1 = new $*1_ltype();
}
%typemap(argout)
  Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>& OUTPUT,
  Eigen::Matrix<DATA_TYPE, Eigen::Dynamic, 1>& OUTPUT,
	Eigen::Matrix<DATA_TYPE, 1, Eigen::Dynamic>& OUTPUT,
  Eigen::Array<DATA_TYPE, Eigen::Dynamic, Eigen::Dynamic>& OUTPUT,
  Eigen::Array<DATA_TYPE, Eigen::Dynamic, 1>& OUTPUT
{
	npy_intp dims[2] = { $1->rows(), $1->cols() };
	PyObject* res = PyArray_New(&PyArray_Type, 2, dims, DATA_TYPECODE, NULL, (void*) $1->data(), 0, NPY_ARRAY_FARRAY, NULL);
	if (!res) { delete $1; SWIG_fail; }
	array_setbase(res, encapsulate($1, clean<$1_ltype>));
	$result = SWIG_Python_AppendOutput($result,res);
}


%enddef    /* %eigen_numpy_typemaps() macro */



/*****************************************************************************
  The following instantiates the typemaps for various underlying data types:
 ****************************************************************************/

%define PREC_INT8_ARRAY    1025 %enddef
%define PREC_INT16_ARRAY   1035 %enddef
%define PREC_INT32_ARRAY   1045 %enddef
%define PREC_INT64_ARRAY   1055 %enddef
%define PREC_INT128_ARRAY  1065 %enddef
%define PREC_FLOAT_ARRAY   1080 %enddef
%define PREC_DOUBLE_ARRAY  1090 %enddef

%eigen_numpy_typemaps(int       , NPY_INT     , PREC_INT32_ARRAY   )
%eigen_numpy_typemaps(long      , NPY_LONG    , PREC_INT64_ARRAY   )
%eigen_numpy_typemaps(float     , NPY_FLOAT   , PREC_FLOAT_ARRAY   )
%eigen_numpy_typemaps(double    , NPY_DOUBLE  , PREC_DOUBLE_ARRAY  )


#endif /* SWIGPYTHON */
