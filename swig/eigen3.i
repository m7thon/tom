/* -*- C++ -*-  (not really, but good for syntax highlighting) */
#ifdef SWIGPYTHON

// numpy.i is from https://github.com/numpy/numpy/tree/master/doc/swig
%include "numpy.i"

%{
#include "Eigen/Core"
%}

/*****************************************************************************
  Some convenience typemaps:
 ****************************************************************************/
%inline %{
	typedef Eigen::Map<Eigen::Matrix<    int, -1, -1>, 0, Eigen::Stride<-1, -1> >  MatrixMi;
	typedef Eigen::Map<Eigen::Matrix<   long, -1, -1>, 0, Eigen::Stride<-1, -1> >  MatrixMl;
	typedef Eigen::Map<Eigen::Matrix<  float, -1, -1>, 0, Eigen::Stride<-1, -1> >  MatrixMf;
	typedef Eigen::Map<Eigen::Matrix< double, -1, -1>, 0, Eigen::Stride<-1, -1> >  MatrixMd;
	typedef Eigen::Map<Eigen::Array<     int, -1, -1>, 0, Eigen::Stride<-1, -1> >  ArrayMMi;
	typedef Eigen::Map<Eigen::Array<    long, -1, -1>, 0, Eigen::Stride<-1, -1> >  ArrayMMl;
	typedef Eigen::Map<Eigen::Array<   float, -1, -1>, 0, Eigen::Stride<-1, -1> >  ArrayMMf;
	typedef Eigen::Map<Eigen::Array<  double, -1, -1>, 0, Eigen::Stride<-1, -1> >  ArrayMMd;
	typedef Eigen::Map<Eigen::Array<     int, -1,  1>, 0, Eigen::Stride< 0, -1> >   ArrayMi;
	typedef Eigen::Map<Eigen::Array<    long, -1,  1>, 0, Eigen::Stride< 0, -1> >   ArrayMl;
	typedef Eigen::Map<Eigen::Array<   float, -1,  1>, 0, Eigen::Stride< 0, -1> >   ArrayMf;
	typedef Eigen::Map<Eigen::Array<  double, -1,  1>, 0, Eigen::Stride< 0, -1> >   ArrayMd;
%}

namespace Eigen {
    const int Dynamic = -1;
    template<typename T, typename D1, typename D2> class Matrix;
    template<typename T, typename D1, typename D2> class Array;
    typedef Matrix<   int, -1, -1> MatrixXi;
    typedef Matrix< float, -1, -1> MatrixXf;
    typedef Matrix<double, -1, -1> MatrixXd;
    typedef Matrix<   int, -1,  1> VectorXi;
    typedef Matrix< float, -1,  1> VectorXf;
    typedef Matrix<double, -1,  1> VectorXd;
    typedef Matrix<   int,  1, -1> RowVectorXi;
    typedef Matrix< float,  1, -1> RowVectorXf;
    typedef Matrix<double,  1, -1> RowVectorXd;
    typedef  Array<   int, -1, -1> ArrayXXi;
    typedef  Array< float, -1, -1> ArrayXXf;
    typedef  Array<double, -1, -1> ArrayXXd;
    typedef  Array<   int, -1,  1> ArrayXi;
    typedef  Array< float, -1,  1> ArrayXf;
typedef  Array<double, -1,  1> ArrayXd;
}

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
  The actual typemaps follow:
 ****************************************************************************/

%define %array_conversion_code(DATA_TYPECODE, ROWS, COLS, DIM)
    PyArrayObject * ary = obj_to_array_no_conversion($input, DATA_TYPECODE);
    if (ary == NULL) SWIG_fail;
    if (array_numdims(ary) != DIM) { PyErr_SetString(PyExc_ValueError, "array must be " #DIM "-dimensional"); SWIG_fail; }
    int rows = array_size(ary, 0);
    int cols = DIM == 1 ? 1 : array_size(ary, 1);
    if (DIM == 1 and ROWS == 1 and COLS != 1) { cols = rows; rows = 1; }
    if (ROWS != -1 and ROWS != rows) { PyErr_SetString(PyExc_ValueError, "array must have exactly" #ROWS "rows"); SWIG_fail; }
    if (COLS != -1 and COLS != cols) { PyErr_SetString(PyExc_ValueError, "array must have exactly" #COLS "columns"); SWIG_fail; }
    int inner = array_stride(ary,0)/PyArray_ITEMSIZE(ary);
    int outer = DIM == 1 ? 0 : array_stride(ary,1)/PyArray_ITEMSIZE(ary);
%enddef

%define %eigen_standard_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, EIGEN_TYPE, ROWS, COLS, DIM)
    %typemap(typecheck, precedence = DATA_TYPEPRECEDENCE, fragment="Eigen_NumPy_Utilities")
      EIGEN_TYPE<DATA_TYPE, ROWS, COLS>,
      const EIGEN_TYPE<DATA_TYPE, ROWS, COLS> &
    {
      $1 = is_array($input) && type_match(array_type($input), DATA_TYPECODE);
    }

    %typemap(in, fragment="Eigen_NumPy_Utilities")
      const EIGEN_TYPE<DATA_TYPE, ROWS, COLS> & (EIGEN_TYPE<DATA_TYPE, ROWS, COLS> temp)
    {
        %array_conversion_code(DATA_TYPECODE, ROWS, COLS, DIM);
        temp = Eigen::Map<Eigen::Matrix<DATA_TYPE, -1, -1>, 0, Eigen::Stride<-1, -1> >
               ((DATA_TYPE*) array_data(ary), rows, cols, Eigen::Stride<-1, -1>(outer, inner));
        $1 = &temp;
    }

    %typemap(in, optimal=1, fragment="Eigen_NumPy_Utilities")
      EIGEN_TYPE<DATA_TYPE, ROWS, COLS>
    {
        %array_conversion_code(DATA_TYPECODE, ROWS, COLS, DIM);
        $1 = Eigen::Map<Eigen::Matrix<DATA_TYPE, -1, -1>, 0, Eigen::Stride<-1, -1> >
             ((DATA_TYPE*) array_data(ary), rows, cols, Eigen::Stride<-1, -1>(outer, inner));
    }

    // (Output) & [does not pass ownership]
    %typemap(out, fragment="Eigen_NumPy_Utilities")
      EIGEN_TYPE<DATA_TYPE, ROWS, COLS> &
    {
        npy_intp dims[DIM]; dims[0] = $1->rows(); if (DIM == 2) { dims[1] = $1->cols(); }
        $result = PyArray_New(&PyArray_Type, DIM, dims, DATA_TYPECODE, NULL, (void*) const_cast<$1_ltype>($1)->data(), 0, NPY_ARRAY_FARRAY, NULL);
        if (!$result) SWIG_fail;
        array_setbase($result,encapsulate($1, NULL));
    }

    // (Output) const & [does not pass ownership]
    %typemap(out, fragment="Eigen_NumPy_Utilities")
      const EIGEN_TYPE<DATA_TYPE, ROWS, COLS> &
    {
        npy_intp dims[DIM]; dims[0] = $1->rows(); if (DIM == 2) { dims[1] = $1->cols(); }
        $result = PyArray_New(&PyArray_Type, DIM, dims, DATA_TYPECODE, NULL, (void*) const_cast<$1_ltype>($1)->data(), 0, NPY_ARRAY_FARRAY_RO, NULL);
        if (!$result) SWIG_fail;
        array_setbase($result,encapsulate($1, NULL));
    }

    // (Output) * [passed ownership]
    %typemap(out, fragment="Eigen_NumPy_Utilities")
      EIGEN_TYPE<DATA_TYPE, ROWS, COLS> *
    {
        npy_intp dims[DIM]; dims[0] = $1->rows(); if (DIM == 2) { dims[1] = $1->cols(); }
        $result = PyArray_New(&PyArray_Type, DIM, dims, DATA_TYPECODE, NULL, (void*) $1->data(), 0, NPY_ARRAY_FARRAY, NULL);
        if (!$result) SWIG_fail;
        array_setbase($result, encapsulate($1, clean<$1_ltype>));
    }

    // (Output) std::shared_ptr
    %typemap(out, optimal = "1", fragment="Eigen_NumPy_Utilities")
      std::shared_ptr<EIGEN_TYPE<DATA_TYPE, ROWS, COLS> >
    {
        auto c_obj = new $1_ltype($1);
        npy_intp dims[DIM]; dims[0] = c_obj->rows(); if (DIM == 2) { dims[1] = c_obj->cols(); }
        $result = PyArray_New(&PyArray_Type, DIM, dims, DATA_TYPECODE, NULL, (void*) c_obj->data(), 0, NPY_ARRAY_FARRAY, NULL);
        if (!$result) { delete c_obj; SWIG_fail; }
        array_setbase($result, encapsulate(c_obj, clean<$1_ltype*>));
    }

    // (Output) [pass by value: moves object and passes ownership]
    %typemap(out, optimal = "1", fragment="Eigen_NumPy_Utilities")
      EIGEN_TYPE<DATA_TYPE, ROWS, COLS>
    {
        auto c_obj = new $1_ltype($1);
        npy_intp dims[DIM]; dims[0] = c_obj->rows(); if (DIM == 2) { dims[1] = c_obj->cols(); }
        $result = PyArray_New(&PyArray_Type, DIM, dims, DATA_TYPECODE, NULL, (void*) c_obj->data(), 0, NPY_ARRAY_FARRAY, NULL);
        if (!$result) { delete c_obj; SWIG_fail; }
        array_setbase($result, encapsulate(c_obj, clean<$1_ltype*>));
    }

    // (Argout) & OUTPUT [creates new object, passes ownership]
    %typemap(in, numinputs=0)
      EIGEN_TYPE<DATA_TYPE, ROWS, COLS> & OUTPUT
    {
        $1 = new $*1_ltype();
    }
    %typemap(argout)
      EIGEN_TYPE<DATA_TYPE, ROWS, COLS> & OUTPUT
    {
        npy_intp dims[DIM]; dims[0] = $1->rows(); if (DIM == 2) { dims[1] = $1->cols(); }
        PyObject* res = PyArray_New(&PyArray_Type, DIM, dims, DATA_TYPECODE, NULL, (void*) $1->data(), 0, NPY_ARRAY_FARRAY, NULL);
        if (!res) { delete $1; SWIG_fail; }
        array_setbase(res, encapsulate($1, clean<$1_ltype>));
        $result = SWIG_Python_AppendOutput($result,res);
    }
%enddef

%define %eigen_map_input_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, EIGEN_BASE, EIGEN_TYPE, ROWS, COLS, DIM, STRIDE)
    %typemap(typecheck, precedence = DATA_TYPEPRECEDENCE, fragment="Eigen_NumPy_Utilities")
    const EIGEN_BASE<Eigen::Map<EIGEN_TYPE<DATA_TYPE, ROWS, COLS>, 0, STRIDE > >&
    {
        $1 = is_array($input) && type_match(array_type($input), DATA_TYPECODE);
    }

    %typemap(in, fragment="Eigen_NumPy_Utilities")
    const EIGEN_BASE<Eigen::Map<EIGEN_TYPE<DATA_TYPE, ROWS, COLS>, 0, STRIDE > >&
    {
        %array_conversion_code(DATA_TYPECODE, ROWS, COLS, DIM);
        $1 = new Eigen::Map<EIGEN_TYPE<DATA_TYPE, ROWS, COLS>, 0, STRIDE > ((DATA_TYPE*) array_data(ary), rows, cols, STRIDE(outer, inner));
    }

    %typemap(freearg)
    const EIGEN_BASE<Eigen::Map<EIGEN_TYPE<DATA_TYPE, ROWS, COLS>, 0, STRIDE > >&
    { if ($1) delete $1; }
%enddef

%define %eigen_other_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, EIGEN_BASE, EIGEN_TYPE, ROWS, COLS, DIM)
    %typemap(typecheck, precedence = DATA_TYPEPRECEDENCE, fragment="Eigen_NumPy_Utilities")
      const EIGEN_BASE<EIGEN_TYPE<DATA_TYPE, ROWS, COLS> >&
    {
      $1 = is_array($input) && type_match(array_type($input), DATA_TYPECODE);
    }

    // (Input) const & Eigen::Derived [makes copy, passes Eigen::Derived object]
    %typemap(in, fragment="Eigen_NumPy_Utilities")
      const EIGEN_BASE<EIGEN_TYPE<DATA_TYPE, ROWS, COLS> >& (EIGEN_TYPE<DATA_TYPE, ROWS, COLS> temp)
    {
        %array_conversion_code(DATA_TYPECODE, ROWS, COLS, DIM);
        temp = Eigen::Map<EIGEN_TYPE<DATA_TYPE, -1, -1>, 0, Eigen::Stride<-1, -1> >
               ((DATA_TYPE*) array_data(ary), rows, cols, Eigen::Stride<-1, -1>(outer, inner));
        $1 = &temp;
    }

    // MARK: (Argout) const & OUTPUT [creates new object, passes ownership]
    %typemap(in, numinputs=0)
      const EIGEN_BASE<EIGEN_TYPE<DATA_TYPE, ROWS, COLS> >& OUTPUT
    {
        $1 = new EIGEN_TYPE<DATA_TYPE, ROWS, COLS>();
    }
    %typemap(argout)
      const EIGEN_BASE<EIGEN_TYPE<DATA_TYPE, ROWS, COLS> >& OUTPUT
    {
        npy_intp dims[DIM]; dims[0] = $1->rows(); if (DIM == 2) { dims[1] = $1->cols(); }
        PyObject* res = PyArray_New(&PyArray_Type, DIM, dims, DATA_TYPECODE, NULL, (void*) $1->derived().data(), 0, NPY_ARRAY_FARRAY, NULL);
        if (!res) { delete $1; SWIG_fail; }
        array_setbase(res, encapsulate($1, clean<EIGEN_TYPE<DATA_TYPE,ROWS,COLS>*>));
        $result = SWIG_Python_AppendOutput($result,res);
    }
%enddef

/*****************************************************************************
  The following instantiates the typemaps for various underlying data types:
 ****************************************************************************/

%define %eigen_templated_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, EIGEN_BASE)
    %eigen_map_input_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, EIGEN_BASE, Eigen::Matrix, -1, -1, 2, %arg(Eigen::Stride<-1, -1>))
    %eigen_map_input_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, EIGEN_BASE, Eigen::Array, -1, -1, 2, %arg(Eigen::Stride<-1, -1>))
    %eigen_map_input_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, EIGEN_BASE, Eigen::Array, -1,  1, 1, %arg(Eigen::Stride< 0, -1>))

    %eigen_other_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, EIGEN_BASE, Eigen::Matrix, -1, -1, 2)
    %eigen_other_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, EIGEN_BASE, Eigen::Matrix, -1,  1, 2)
    %eigen_other_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, EIGEN_BASE, Eigen::Matrix,  1, -1, 2)
    %eigen_other_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, EIGEN_BASE,  Eigen::Array, -1, -1, 2)
    %eigen_other_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, EIGEN_BASE,  Eigen::Array, -1,  1, 1)
%enddef

%define %eigen_numpy_typemaps(DATA_TYPE, DATA_TYPECODE, DATA_TYPEPRECEDENCE)
    %eigen_standard_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, Eigen::Matrix, -1, -1, 2)
    %eigen_standard_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, Eigen::Matrix, -1,  1, 2)
    %eigen_standard_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, Eigen::Matrix,  1, -1, 2)
    %eigen_standard_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE,  Eigen::Array, -1, -1, 2)
    %eigen_standard_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE,  Eigen::Array, -1,  1, 1)

    %eigen_templated_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, Eigen::MatrixBase)
    %eigen_templated_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, Eigen::ArrayBase)
    %eigen_templated_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, Eigen::DenseBase)
    %eigen_templated_typemaps(DATA_TYPE, DATA_TYPEPRECEDENCE, DATA_TYPECODE, Eigen::EigenBase)
%enddef

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
