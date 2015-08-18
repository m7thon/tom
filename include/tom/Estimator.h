#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include "tom.h"

namespace tom {

double ApproxNormalCDFInverse(double p) {
	constexpr double c[] = {2.515517, 0.802853, 0.010328};
	constexpr double d[] = {1.432788, 0.189269, 0.001308};
	if (p < 0.5) {
		double t = std::sqrt(-2.0*std::log(p));
		return -( t - ((c[2]*t + c[1])*t + c[0]) / (((d[2]*t + d[1])*t + d[0])*t + 1.0) );
	}
	else {
		double t = std::sqrt(-2.0*std::log(1-p));
		return t - ((c[2]*t + c[1])*t + c[0]) / (((d[2]*t + d[1])*t + d[0])*t + 1.0);
	}
}


/**
 * This class computes estimates \f$\hat{f}(\bar{x})\f$ and corresponding variance estimates for sequences \f$\bar{x}\f$ based on a suffix tree representation of a sample sequence.
*/
class Estimator {
    public:
	///** create an uninitialized \c Estimator. */
	//Estimator() : pos_() {};

	/** create an \c Estimator from a given \c sfxTree -- a suffix tree representation of a sample sequence */
	Estimator(const std::shared_ptr<stree::STree>& sfxTree) : s_(sfxTree), stree_(sfxTree) {
		s_.pos_ = stree::Position(stree_);
		nO_ = stree_->sequence().nO(); nU_ = stree_->sequence().nU();
		len_ = ( stree_->sequence().length() );
		uProbs_ = Eigen::VectorXd::Ones(std::max(1, nU_));
		for (Symbol u = 0; u < nU_; ++u) {
			s_.pos_ = stree::Position(stree_); s_.pos_.toSymbol(u);
			uProbs_(u) = (double)(s_.pos_.count()) / len_;
		}
		s_.pos_ = stree::Position(stree_);
	}

	/** return the size of the underlying input alphabet */
	unsigned int nU() { return nU_; }

	/** return the size of the underlying output alphabet */
	unsigned int nO() { return nO_; }

  /** return the current estimate */
	double f(const Sequence& seq) { estimateVariance_ = false; reset(); extendBy(seq); eval(); return s_.f_; }

	/** return a variance estimate for the current estiamte */
	double v(const Sequence& seq) { estimateVariance_ = true; reset(); extendBy(seq); eval(); return s_.v_; }

  SWIGCODE(%apply Eigen::MatrixXd& OUTPUT { Eigen::MatrixXd& F };)
  SWIGCODE(%apply Eigen::MatrixXd& OUTPUT { Eigen::MatrixXd& V };)
  /** return (in the output argument \c F) the matrix of estimates \f$\hat{F}^{I,J}=[\hat{f}(\bar{x}_j\bar{x}_i)]_{i,j}\f$ for the given set \c chaSeqs of characteristic sequences \f$\bar{x}_i\f$ and the set \c indSeqs of indicative sequences \f$\bar{x}_j\f$. */
	void f(Eigen::MatrixXd& F, const Sequences& chaSeqs, const Sequences& indSeqs) {
		unsigned int rows = chaSeqs.size(), cols = indSeqs.size();
		F.resize(rows, cols);
		estimateVariance_ = false;
		for (unsigned int j = 0; j < cols; ++j) {
			reset(); extendBy(indSeqs[j]);
			State s = s_;
			for (unsigned int i = 0; i < rows; ++i) {
				reset(s); extendBy(chaSeqs[i]); eval();
				F.coeffRef(i,j) = s_.f_;
			}
		}
	}

	/** return (in the output argument \c F) the matrix of estimates \f$\hat{F}_z^{I,J}=[\hat{f}(\bar{x}_jz\bar{x}_i)]_{i,j}\f$ for the given set \c chaSeqs of characteristic sequences \f$\bar{x}_i\f$, the set \c indSeqs of indicative sequences \f$\bar{x}_j\f$ and the input-output pair \c z = (\c u,\c o). */
	void f(Eigen::MatrixXd& F, const Sequences& chaSeqs, const Sequences& indSeqs, Symbol o, Symbol u = 0) {
		unsigned int rows = chaSeqs.size(), cols = indSeqs.size();
		F.resize(rows, cols);
		estimateVariance_ = false;
		for (unsigned int j = 0; j < cols; ++j) {
			reset(); extendBy(indSeqs[j]); extendBy(o,u);
			State s = s_;
			for (unsigned int i = 0; i < rows; ++i) {
				reset(s); extendBy(chaSeqs[i]); eval();
				F.coeffRef(i,j) = s_.f_;
			}
		}
	}

	/** return (in the output argument \c V) the matrix of variance estimates \f$\hat{V}^{I,J}[\hat{\rm{Var}}(\hat{f}(\bar{x}_j\bar{x}_i))]_{i,j}\f$ corresponding to the estimates \f$\hat{f}(\bar{x}_j\bar{x}_i)\f$ for the given set \c chaSeqs of characteristic sequences \f$\bar{x}_i\f$ and the set \c indSeqs of indicative sequences \f$\bar{x}_j\f$. */
	void v(Eigen::MatrixXd& V, const Sequences& chaSeqs, const Sequences& indSeqs) {
		unsigned int rows = chaSeqs.size(), cols = indSeqs.size();
		V.resize(rows, cols);
		estimateVariance_ = true;
		for (unsigned int j = 0; j < cols; ++j) {
			reset(); extendBy(indSeqs[j]);
			State s = s_;
			for (unsigned int i = 0; i < rows; ++i) {
				reset(s); extendBy(chaSeqs[i]); eval();
				V.coeffRef(i,j) = s_.v_;
			}
		}
	}

	/** return (in the output argument \c V) the matrix of variance estimates \f$\hat{V}^{I,J}[\hat{\rm{Var}}(\hat{f}(\bar{x}_jz\bar{x}_i))]_{i,j}\f$ corresponding to the estimates \f$\hat{f}(\bar{x}_jz\bar{x}_i)\f$ for the given set \c chaSeqs of characteristic sequences \f$\bar{x}_i\f$, the set \c indSeqs of indicative sequences \f$\bar{x}_j\f$, and the input-output pair \c z = (\c u,\c o). */
	void v(Eigen::MatrixXd& V, const Sequences& chaSeqs, const Sequences& indSeqs, Symbol o, Symbol u = 0) {
		unsigned int rows = chaSeqs.size(), cols = indSeqs.size();
		V.resize(rows, cols);
		estimateVariance_ = true;
		for (unsigned int j = 0; j < cols; ++j) {
			reset(); extendBy(indSeqs[j]); extendBy(o,u);
			State s = s_;
			for (unsigned int i = 0; i < rows; ++i) {
				reset(s); extendBy(chaSeqs[i]); eval();
				V.coeffRef(i,j) = s_.v_;
			}
		}
	}

	/** return (in the output argument \c F) the matrix of estimates \f$\hat{F}_z^{I,J}=[\hat{f}(\bar{x}_j\bar{x}_i)]_{i,j}\f$, and (in the output argument \c V) the corresponding matrix of variance estimates \f$\hat{V}^{I,J}[\hat{\rm{Var}}(\hat{f}(\bar{x}_j\bar{x}_i))]_{i,j}\f$ for the given set \c chaSeqs of characteristic sequences \f$\bar{x}_i\f$ and the set \c indSeqs of indicative sequences \f$\bar{x}_j\f$. */
	void fv(Eigen::MatrixXd& F, Eigen::MatrixXd& V, const Sequences& chaSeqs, const Sequences& indSeqs) {
		unsigned int rows = chaSeqs.size(), cols = indSeqs.size();
		F.resize(rows, cols); V.resize(rows, cols);
		estimateVariance_ = true;
		for (unsigned int j = 0; j < cols; ++j) {
			reset(); extendBy(indSeqs[j]);
			State s = s_;
			for (unsigned int i = 0; i < rows; ++i) {
				reset(s); extendBy(chaSeqs[i]); eval();
				F.coeffRef(i,j) = s_.f_;
				V.coeffRef(i,j) = s_.v_;
			}
		}
	}

	/** return (in the output argument \c F) the matrix of estimates \f$\hat{F}_z^{I,J}=[\hat{f}(\bar{x}_jz\bar{x}_i)]_{i,j}\f$, and (in the output argument \c V) the corresponding matrix of variance estimates \f$\hat{V}^{I,J}[\hat{\rm{Var}}(\hat{f}(\bar{x}_jz\bar{x}_i))]_{i,j}\f$ for the given set \c chaSeqs of characteristic sequences \f$\bar{x}_i\f$, the set \c indSeqs of indicative sequences \f$\bar{x}_j\f$ and the input-output pair \c z = (\c u,\c o). */
	void fv(Eigen::MatrixXd& F, Eigen::MatrixXd& V, const Sequences& chaSeqs, const Sequences& indSeqs, Symbol o, Symbol u = 0) {
		unsigned int rows = chaSeqs.size(), cols = indSeqs.size();
		F.resize(rows, cols); V.resize(rows, cols);
		estimateVariance_ = true;
		for (unsigned int j = 0; j < cols; ++j) {
			reset(); extendBy(indSeqs[j]); extendBy(o,u);
			State s = s_;
			for (unsigned int i = 0; i < rows; ++i) {
				reset(s); extendBy(chaSeqs[i]); eval();
				F.coeffRef(i,j) = s_.f_;
				V.coeffRef(i,j) = s_.v_;
			}
		}
	}

  SWIGCODE(%clear Eigen::MatrixXd& F;)
  SWIGCODE(%clear Eigen::MatrixXd& V;)

	int nO_;              ///< the size of the output alphabet
	int nU_;              ///< the size of the input alphabet
	unsigned int len_;    ///< the length of the sample sequence from which the estimates are calculated
	Eigen::VectorXd uProbs_; ///< the input Symbol probabilities for the case of an iid ("blind") input policy.\ These are estimated on constructing the \c Estimator and may be overwritten if they are known exactly.

	double nPseudoCounts_ = 1;
	double zConfidenceIntervalSize_ = 0;
	double addToVariance_ = 0;
	double minimumVariance_ = 0;
	double applyExponentToVariance_ = 1;

private:
  struct State {
      State(const std::shared_ptr<stree::STree> stree) : pos_(stree) {}
		stree::Position pos_; ///< the position in the suffix tree for the currently estimated sequence.
		double f_  = 1;       ///< related to the current estimate (used internally in different ways)
		double v_  = 1;
		double fB_ = 1;       ///< related to the current Bayesian estimate (used internally in different ways)
	};


  State s_;
	bool estimateVariance_ = true;

	const std::shared_ptr<stree::STree> stree_; ///< a pointer to the underlying \c STree

	void reset() { s_.pos_ = stree::Position(stree_); s_.f_ = s_.fB_ = s_.v_ = 1; }
	void reset(const State& s) { s_ = s; }

	void extendBy(Symbol o, Symbol u = 0) {
		if (nU_ == 0) s_.pos_.toSymbol(o);
		else {
			s_.pos_.toSymbol(u); double cu = s_.pos_.count();
			s_.pos_.toSymbol(o); double co = s_.pos_.count();
			if (cu == 0) s_.f_ *= ( (double)(1) / nO_ );
			else         s_.f_ *= ( (double)(co) / cu );
			if (estimateVariance_) {
				if (cu < 2) { s_.v_ *= 0.75; }
				else {
					double pB = ( co + nPseudoCounts_ ) / ( cu + nPseudoCounts_ );
					s_.fB_ *= pB;
					double pBv = ( co + 0.5 * nPseudoCounts_ ) / ( cu + nPseudoCounts_ );
					double vB = pBv * ( 1 - pBv ) / ( cu );
					s_.v_ *= std::max( pB * pB - vB, (double)(0) );
				}
			}
		}
	}

	void extendBy(const Sequence& seq) {
		if (nU_ == 0) s_.pos_.toSequence(seq);
		else for (unsigned int i = 0; i < seq.length(); ++i) extendBy(seq.o(i), seq.u(i));
	}

	void eval() {
		if (nU_ == 0) {
			double n = len_ - s_.pos_.depth() + 1;
			double x = s_.pos_.count();
			s_.f_ = x / n;
			if (estimateVariance_) {
				if (n < 2) { s_.v_ = 0.25; }
				else {
					// Compute a confidence interval for the estimate f, and use the bound that is closer to 1/2 to estimate the variance:
					double fB = ( x + 0.5 * nPseudoCounts_ ) / ( n + nPseudoCounts_ );
					double confidenceIntervalRadius = zConfidenceIntervalSize_ * std::sqrt( fB * ( 1 - fB ) / ( n + nPseudoCounts_ ) );
					if ( fB <= 0.5 ) fB = std::min( fB + confidenceIntervalRadius, (double)(0.5) );
					else fB = std::max( (double)(0.5), fB - confidenceIntervalRadius );
					s_.v_ = fB * ( 1 - fB ) / ( n );
				}
				s_.v_ = std::pow( std::max( s_.v_ + addToVariance_, minimumVariance_ ), applyExponentToVariance_ );
			}
		// 		else if (useExactCI_ == true){
		// 			// Compute a Clopper-Pearson "exact" confidence interval for the estimate, and use the bound closer to 1/2 to estimate the variance:
		// 			double fB = s_.f_;
		// 			if (fB < 0.5) {
		// 				double a = x+1; double b = n-x;
		// 				double log_beta = std::lgamma(a) + std::lgamma(b) - std::lgamma(a+b);
		// 				int error = 0;
		// 				fB = xinbta(a, b, log_beta, 1-maxProbAlphaError_, error);
		// 				if (fB > 0.5) fB = 0.5;
		// 			}
		// 			else if (fB > 0.5) {
		// 				double a = x; double b = n-x+1;
		// 				double log_beta = std::lgamma(a) + std::lgamma(b) - std::lgamma(a+b);
		// 				int error = 0;
		// 				fB = xinbta(a, b, log_beta, maxProbAlphaError_, error);
		// 				if (fB < 0.5) fB = 0.5;
		// 			}
		// 			s_.v_ = std::max( fB * ( 1 - fB ) / (n), (double)(1) / ( (double)(n) * (double)(n) ) );
		// 		}
		// 	}
		}
		else if (estimateVariance_) { // nU != 0
			s_.v_ = std::pow( (s_.v_ == 1 ? (double)(1) / (double)(len_) * (double)(len_) : s_.fB_ * s_.fB_ - s_.v_), applyExponentToVariance_ );
		}
	}

}; // class Estimator

#ifdef SWIG
%pythoncode %{
	EstimatorVarianceParams = {
		'nPseudoCounts_': 1, 'zConfidenceIntervalSize_': 1, 'addToVariance_': 0,
		'minimumVariance_': 1e-15, 'applyExponentToVariance_': 0.5 }
	EstimatorExactVarianceParams = {
		'nPseudoCounts_': 0, 'zConfidenceIntervalSize_': 0, 'addToVariance_': 0,
		'minimumVariance_': 0, 'applyExponentToVariance_': 1 }
%}

#endif

} // namespace tom

#endif // ESTIMATOR_H
