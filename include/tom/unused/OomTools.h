#ifndef OOM_TOOLS_H_
#define OOM_TOOLS_H_

#include "../tom.h"

namespace tom {

SHARED_PTR<Oom> hmmToOom(const Eigen::MatrixXd& T, const Eigen::MatrixXd& E, const Eigen::VectorXd& w, bool transition_first = false); 


#ifndef SWIG
// ##################################################################
//                         IMPLEMENTATION
// ##################################################################
SHARED_PTR<Oom> hmmToOom(const Eigen::MatrixXd& T, const Eigen::MatrixXd& E, const Eigen::VectorXd& w, bool transition_first) {
	SHARED_PTR<Oom> oom(new Oom());
	oom.setSize(T.rows(), E.cols(), 0);
	nU_ = 0;
	nO_ = E.cols();
	dim_ = T.rows();
	oom.sig() = Eigen::VectorXd::Ones(oom.dim());
	for (Symbol o = 0; o < oom.nO(); ++o) {
		if (transition_first) oom.tau(o,0) = E.col(o).asDiagonal() * T.transpose(); 
		else oom.tau(o,0) = T.transpose() * E.col(o).asDiagonal(); 
	}
	oom.w0() = w;
	oom.init();
	return oom;
}

#endif // SWIG

} // namespace tom

#endif // OOM_TOOLS_H_
