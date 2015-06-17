/* -*- C -*-  (not really, but good for syntax highlighting) */
%module eigen3
%naturalvar;

%{ 
#define SWIG_FILE_WITH_INIT
#include <iostream>
#include "../include/external/Eigen/Core"
	using namespace std;
	
	double sum(const Eigen::Ref<Eigen::MatrixXd>& m) {
		return m.sum();
	}
	template<typename T>
	double sum2(const Eigen::MatrixBase<T>&m ) {
		return m.sum();
	}
%} 

%include "eigen3.i"

%init %{
  import_array();
%}

double sum(const Eigen::Ref<Eigen::MatrixXd>& m);
template<typename T>
double sum2(const Eigen::MatrixBase<T>&m );
%template(sum2) sum2<MatrixMd>;

//%include "../Macros.h"
//TEMPLATE(C,Learner::C,MatrixMd)


