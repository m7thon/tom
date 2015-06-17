/**
 * @file   OomTools.h
 * @author Michael Thon
 * 
 * @brief  This file provides tools for working with OOMs.
 */

#ifndef OOM_TOOLS_H_
#define OOM_TOOLS_H_

namespace tom {

/**
 * Create an OOM from a given HMM. The HMM parameters must be specified as follows:
 * @param T the state transition probabilities: \f$T_{i,j} = P(s_{t+1}=j | s_{t}=i)\f$
 * @param E the emission probabilities: \f$E_{i,j} = P(o_t=j | s_t=i)\f$
 * @param w the initial state probabilities: \f$w_i = P(s_0 = i)\f$.
 */
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
