/**
 * @file   Learner.h
 * @author Michael Thon <mthon@jacobs-university.de>
 * @date   Mon Aug  2 14:30:36 2010
 * 
 * @brief  Functions for learning OOMs from data
 */

#ifndef _OOM_TRAIN_H_
#define _OOM_TRAIN_H_

#include <cstdlib>

namespace tom {

class Learner {
public:
	Learner() {};
	enum { NO_WEIGHTS = 0, RC_WEIGHTS = 1, FULL_WEIGHTS = 2};
	enum { EC_ALGO, SVD_ALGO };
	enum { NO_CACHE = 0, CACHE_F = 1, CACHE_W = 2, CACHE_FZ = 4, CACHE_WZ = 8, CACHE_CQ = 16 };
	int nU_;       /**< the size of the input alphabet */
  int nO_;       /**< the size of the output alphabet */
	int dim_;      /**< the target dimension of the \a Oom model to be learnt */
	unsigned long seqLen_; /**< the length of the sample sequence to use for learning */
	int weightingScheme_; /* how to use weights is the learning procedure -- if permitted by the chosen algorithm. The possible choices are:
													 - NO_WEIGHTS: Don't use weights at all
													 - FULL_WEIGHTS: Use all available weights
													 - RC_WEIGHTS: Use only estimated row and column weights
													 - RCREDUCED_WEIGHTS: Use row and column weights that are computed from the full weights */
	int algorithm_; /**< which algorithm to use to estimate the \a Oom model. Currently, the only choices are:
										    - SVD_ALGO: Spectral learning
												- EC_ALGO: Same as spectral learning, but the SVD is computed via EM */
	int cache_;     /**< determines which values in the computation to cache */
	double ec_err_threshold_;
	int ec_max_iterations_;
	int ec_iterations_;
	double ec_err_;

	Sequence trainSequence_;    /**< the sample sequence from which to learn the \a Oom model. */
	stree::STree suffixTree_;   /**< a suffix tree representation of the \a trainSequence_ or a subsequence */ 
	Estimator estimator_;       /**< the \a Estimator to use to obtain estimates from the sample sequence. */
	
	SHARED_PTR<Sequences> characteristicSequences_;   /**< the set of characteristic sequences */
	SHARED_PTR<Sequences> indicativeSequences_;       /**< the set of indicative sequences */
		
/** @name Accessors */
//@{
	SHARED_PTR<Sequences> characteristicSequences() { return characteristicSequences_; }
	void characteristicSequences(SHARED_PTR<Sequences> characteristicSequences_new) { characteristicSequences_ = characteristicSequences_new; }

	SHARED_PTR<Sequences> indicativeSequences() { return indicativeSequences_; }
	void indicativeSequences(SHARED_PTR<Sequences> indicativeSequences_new) { indicativeSequences_ = indicativeSequences_new; }
	
	Eigen::MatrixXd& C() { return C_; }
	template< typename D >
	void C(const Eigen::MatrixBase<D>& C_new) { C_ = C_new; }
	Eigen::MatrixXd& Q() { return Q_; }
	template< typename D >
	void Q(const Eigen::MatrixBase<D>& Q_new) { Q_ = Q_new; }

	Eigen::MatrixXd& F() { return F_; }
	template< typename D >
	void F(const Eigen::MatrixBase<D>& F_new) { F_ = F_new; }

	Eigen::MatrixXd& Fz(Symbol o, Symbol u = 0) { return Fz_(o,u); }
	template< typename D >
	void Fz(Symbol o, const Eigen::MatrixBase<D>& Fz_new) { Fz_(o,0) = Fz_new; }
	template< typename D >
	void Fz(Symbol o, Symbol u, const Eigen::MatrixBase<D>& Fz_new) { Fz_(o,u) = Fz_new; }

	Eigen::MatrixXd& FI() { return FI_; }
	template< typename D >
	void FI(const Eigen::MatrixBase<D>& FI_new) { FI_ = FI_new; }

	Eigen::MatrixXd& FJ() { return FJ_; }
	template< typename D >
	void FJ(const Eigen::MatrixBase<D>& FJ_new) { FJ_ = FJ_new; }

	Eigen::MatrixXd& W() { return W_; }
	template< typename D >
	void W(const Eigen::MatrixBase<D>& W_new) { W_ = W_new; }

	Eigen::MatrixXd& WI() { return WI_; }
	template< typename D >
	void WI(const Eigen::MatrixBase<D>& WI_new) { WI_ = WI_new; }

	Eigen::MatrixXd& WJ() { return WJ_; }
	template< typename D >
	void WJ(const Eigen::MatrixBase<D>& WJ_new) { WJ_ = WJ_new; }

	Eigen::MatrixXd& Wz(Symbol o, Symbol u = 0) { return Wz_(o,u); }
	template< typename D >
	void Wz(Symbol o, const Eigen::MatrixBase<D>& Wz_new) { Wz_(o,0) = Wz_new; }
	template< typename D >
	void Wz(Symbol o, Symbol u, const Eigen::MatrixBase<D>& Wz_new) { Wz_(o,u) = Wz_new; }

//@}

/** @name Basic OOM functionality */
//@{
	Oom* weightedSpectral();
	void init();
	void computeCQ();
	void clearCQ() { clear(C_); clear(Q_); }
	Oom* oom();

//@}


	~Learner();
	
private:
	Eigen::MatrixXd C_, Q_;
	Eigen::MatrixXd F_, FI_, FJ_, W_, WI_, WJ_;
	Eigen::Array<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic> Fz_, Wz_;
	Eigen::VectorXd s_;
	bool initialized_;
	void clear(Eigen::MatrixXd& mat) { mat = Eigen::MatrixXd::Zero(0,0); }
	bool know(Eigen::MatrixXd& mat) { return (mat.diagonalSize() != 0); }
	Learner(const Learner&);
	Learner& operator=(const Learner&);
};

TEMPLATE(C,Learner::C,MatrixMd)
TEMPLATE(Q,Learner::Q,MatrixMd)
TEMPLATE(F,Learner::F,MatrixMd)
TEMPLATE(FI,Learner::FI,MatrixMd)
TEMPLATE(FJ,Learner::FJ,MatrixMd)
TEMPLATE(Fz,Learner::Fz,MatrixMd)
TEMPLATE(W,Learner::W,MatrixMd)
TEMPLATE(WI,Learner::WI,MatrixMd)
TEMPLATE(WJ,Learner::WJ,MatrixMd)
TEMPLATE(Wz,Learner::Wz,MatrixMd)



Learner::~Learner() {
}

void Learner::init() {
	Fz_.resize(nO_, std::max(1, nU_));
	Wz_.resize(nO_, std::max(1, nU_));
}

void Learner::computeCQ() {
  std::cerr << "Computing CQ" << std::endl;
	if (!know(F_)) estimator_.f(F_, *characteristicSequences_, *indicativeSequences_);
	if ((weightingScheme_ & FULL_WEIGHTS) and !know(W_)) {
		estimator_.v(W_, *characteristicSequences_, *indicativeSequences_);
	}
	if (algorithm_ == EC_ALGO) {
		PG_INFO("EC");
		int ec_iterations_ = 1;
		double ec_err_ = 1 + ec_err_threshold_;
		Eigen::MatrixXd C_old;
		C_ = Eigen::MatrixXd::Random(dim_, F_.rows());
		Q_ = pinv(C_ * F_);
		while ((ec_err_ > ec_err_threshold_) and (ec_iterations_ < ec_max_iterations_)) {
			PG_MARK(0,0);
			C_old = C_;
			C_ = pinv(F_ * Q_);
			Q_ = pinv(C_ * F_);
			ec_iterations_ += 1;
			ec_err_ = (C_old - C_).array().abs().sum();
		}
		PG_DONE;
		if (!(cache_ & CACHE_F)) clear(F_); 
		return;
	}
	if (algorithm_ == SVD_ALGO) {
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(F_, Eigen::ComputeThinU | Eigen::ComputeThinV);
		C_ = svd.matrixU().leftCols(dim_).transpose();
		s_ = svd.singularValues().head(dim_);
		double tolerance = std::max(F_.rows(), F_.cols()) * s_(0) * Eigen::NumTraits<double>::epsilon();
		for ( long i = 0; i < dim_; ++i) s_.coeffRef(i) = ( s_.coeff(i) > tolerance ? 1.0 / s_.coeff(i) : 0 );
		Q_ = svd.matrixV().leftCols(dim_) * s_.asDiagonal();
		if (!(cache_ & CACHE_F)) clear(F_); 
		return;
	}
	else std::cout << "Cannot compute CQ" << std::endl;
	return;
}

Oom* Learner::weightedSpectral() {
	const bool isIO = estimator_.isIO_;
	Oom* oom = new Oom();
	oom->setSize(dim_, nO_, nU_);
	estimator_.fv(F_, W_, *characteristicSequences_, *indicativeSequences_);
	W_ = W_.array().inverse(); 
	Eigen::MatrixXd B, A, WA, Fz, Wz, Az, WAz;
	weightedLowRankApproximation(B, A, F_, W_, dim_);
	Az.resize(dim_, F_.cols());
	WA.resize(dim_, F_.cols());
	WAz.resize(dim_, F_.cols());
	for (int j = 0; j < F_.cols(); ++j) {
		Eigen::MatrixXd T = (pinv(W_.col(j).cwiseSqrt().asDiagonal() * B) * W_.col(j).cwiseSqrt().asDiagonal()).array().square();
		for (int i = 0; i < dim_; ++i) {
			WA(i,j) = 1. / (T.row(i) * W_.col(j).array().inverse().matrix());
		}
	}
	for (Symbol u = 0; u < std::max(nU_, 1); ++u) {
		for (Symbol o = 0; o < nO_; ++o) {
			estimator_.fv(Fz, Wz, *characteristicSequences_, *indicativeSequences_);
			Wz = Wz.array().inverse();
			for (int j = 0; j < F_.cols(); ++j) {
				Eigen::MatrixXd T = (pinv(Wz.col(j).cwiseSqrt().asDiagonal() * B) * Wz.col(j).cwiseSqrt().asDiagonal());
				Az.col(j) = T * Fz.col(j);
				for (int i = 0; i < dim_; ++i) {
					WAz(i,j) = 1. / (T.row(i).array().square().matrix() * Wz.col(j).array().inverse().matrix());
				}
			}
			oom->tau(u,o) = EW_TLS(A, Az, WA, WAz);
		}
	}
	estimator_.fv(FI_, WI_, *characteristicSequences_, Sequences(1));
	WI_ = WI_.array().inverse();
	oom->w0() = pinv(WI_.cwiseSqrt().asDiagonal() * B) * WI_.cwiseSqrt().asDiagonal() * FI_;
	
	estimator_.fv(FJ_, WJ_, Sequences(1), *indicativeSequences_);
	WJ_ = WJ_.array().inverse();
	oom->sig() = EW_TLS(A, FJ_, WA, WJ_);
	oom->init();
	return oom;
}


Oom* Learner::oom() {
	const bool isIO = estimator_.isIO_;
	
	if (!know(C_) or !know(Q_)) computeCQ();
	if (!know(FI_)) estimator_.f(FI_, *characteristicSequences_, Sequences(1));
	if (!know(FJ_)) estimator_.f(FJ_, Sequences(1), *indicativeSequences_);

	Oom* oom = new Oom();
	oom->setSize(dim_, estimator_.nO_, estimator_.nU_);
	oom->sig() = FJ_ * Q_;
	oom->w0() = C_ * FI_;
	oom->w0() /= (oom->sig() * oom->w0());
	for (Symbol o = 0; o < estimator_.nO_; ++o) {
		for (Symbol u = 0; u < std::max(1, estimator_.nU_); ++u) {
			if (!know(Fz_(o,u))) estimator_.f(Fz_(o,u), *characteristicSequences_, *indicativeSequences_, o, u);
			oom->tau(o,u) = C_ * Fz_(o,u) * Q_;
		}
	}
	oom->init();
	return oom;
}


// template <typename EstimatorT, typename IndSeqsT, typename ChaSeqsT>
// Oom* learnOom(EstimatorT& est, const IndSeqsT& indSeqs, const ChaSeqsT& chaSeqs, int dim, int nO, int nU = 0, bool weighted = true, bool useEM = true) {	
// 	unsigned int Frows = chaSeqs.size();
// 	unsigned int Fcols = indSeqs.size();

// 	Eigen::MatrixXd* F = estimateF(est, indSeqs, chaSeqs);

// 	Eigen::VectorXd     w_I = Eigen::VectorXd::Ones(Frows);
// 	Eigen::RowVectorXd  w_J = Eigen::RowVectorXd::Ones(Fcols);

// 	// Calculate the row and column weightings if desired
// 	if (weighted) {
// 		Eigen::MatrixXd* F_I = estimateF(est, CoreSeqs(), chaSeqs);
// 		Eigen::MatrixXd* F_J = estimateF(est, indSeqs, CoreSeqs());

// //		w_I = (F->rowwise().sum().array() * (F->sum() - F->rowwise().sum().array())).sqrt();
// //		w_J = (F->colwise().sum().array() * (F->sum() - F->colwise().sum().array())).sqrt();
// 		w_I = F_I->array().sqrt();
// 		for (unsigned int i=0; i<Frows; i++)
// 			if (w_I(i) != 0) w_I(i) = 1.0/w_I(i);

// 		w_J = F_J->array().sqrt();	
// 		for (unsigned int j=0; j<Fcols; j++)
// 			if (w_J(j) != 0) w_J(j) = 1.0/w_J(j);

// 		*F = w_I.asDiagonal() * (*F);
// 		*F = (*F) * w_J.asDiagonal();
		
// 		delete F_I;
// 		delete F_J;
// 	}
	
// 	// Obtain the C and Q matrices by the (weighted) spectral learning procedure:
// 	CQPair* cqp;
// 	if (useEM) cqp = compute_CQ_by_SVD_EM(*F, dim);
// 	else       cqp = compute_CQ_by_SVD(*F, dim);
// 	if (weighted) {
// 		cqp->C = cqp->C * w_I.asDiagonal();
// 		cqp->Q =  w_J.asDiagonal() * cqp->Q;
// 	}
// 	delete F;

// 	// Compute the OOM via the learning equation
// 	Oom* oom = learnOom(est, indSeqs, chaSeqs, *cqp, nO, nU);
	
// 	delete cqp;
// 	return oom;
// }
// #ifdef SWIG
// %template(learnOom) learnOom<Estimator, CoreSeqs, CoreSeqs>;
// #endif // SWIG



// Oom* learnOom(const Sequence<Symbol>& seq, int dim, int eventLength, bool weighted = true, bool useEM = true) {
// 	bool IO = (seq.nU != 0);
// 	STree sfxTree(seq, 0, false, (IO ? 2 : 1));
// 	CoreSeqs coreSeqs(&sfxTree, eventLength);
// 	Estimator est(&sfxTree, IO, seq.nO);
// 	return learnOom(est, coreSeqs, coreSeqs, dim, seq.nO, seq.nU, weighted, useEM);
// }


} // namespace tom

#endif // _OOM_TRAIN_H_
