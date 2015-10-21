/**
 * @file   Hmm.h
 * @author Michael Thon
 *
 * @brief  This file provides EM learning for HMMs and POMDPs.
 */

#ifndef HMM_H
#define HMM_H

#include "tom.h"

namespace tom {

SWIGCODE(%feature("director") EMStopCondition;)
class EMStopCondition {
public:
	int maxEMIterations_;
	double minRelativeImprovement_;
	double previousLog2Likelihood_ = std::numeric_limits<double>::infinity();
	EMStopCondition(int maxEMIterations = 100, double minRelativeImprovement = 0.0001) :
		maxEMIterations_(maxEMIterations), minRelativeImprovement_(minRelativeImprovement) {}
	virtual bool operator()(int iteration, double log2Likelihood) {
		return ((1 - log2Likelihood / previousLog2Likelihood_ < minRelativeImprovement_) or (iteration >= maxEMIterations_));
	}
	virtual ~EMStopCondition() {};
};


SWIGCODE(%feature("python:slot", "tp_repr", functype="reprfunc") Hmm::repr;)
/**
 * This class provides provides a rudimentry structure for HMMs and POMDPs. It purpose is to
 * - create randomly initialized HMMs or POMDPs
 * - learn HMMs / POMDPs from data using EM (Baum-Welch)
 * Further operations are available after conversion into an \c Oom.
 */
class Hmm {
	friend class cereal::access;
	friend class Oom;
public:
	Hmm(int nStates, int nObservations, int nInputs = 0, double exponent = 1, const Random& rnd = Random())
		{ setSize(nStates, nObservations, nInputs); if (exponent != 0) randomize(exponent, rnd); else init(); }

	void randomize(double exponent = 1, const Random& rnd = Random());
	bool normalize();
	void init();

	double trainEM(const Sequence& trainSequence, const EMStopCondition& stopCondition = EMStopCondition());

	/** @name Accessors */ //@{
	int nStates()       const { return dim_; }
	int nObservations() const { return nO_; }
	int nInputs()       const { return nU_; }
 	void setSize(int nStates, int nObservations, int nInputs, bool zeroParameters = false);
	VectorXd& pi()                             { return pi_;      }
	void pi(const VectorXd& _pi)               { pi_ = _pi;       }
	MatrixXd& T(int a = 0 )                    { return T_(a);    }
	void T(const MatrixXd& _T)                 { T_(0) = _T;      }
	void T(int a, const MatrixXd& _Ta)         { T_(a) = _Ta;     }
	VectorXd& E(int o, int a = 0)              { return E_(o,a);  }
	void E(int o, int a, const VectorXd& _Eoa) { E_(o, a) = _Eoa; }
	void E(int o, const VectorXd& _Eo)         { E_(o,0) = _Eo;   }

	MatrixXd& Theta(int o, int a = 0 )         { return Theta_(o,a); }
  //@}

	INSERT_JSON_IO_FUNCTIONS()
	/** return a representation to display in interactive python. */
	std::string repr() const { std::stringstream os; os << "HMM" << " nU: " << nU_ << " nO: " << nO_ << " dim: " << dim_; return os.str(); }

private:
	int dim_, nO_, nU_, nU1_;
	VectorXd pi_;
	Array<MatrixXd, 1, Dynamic> T_;
	Array<VectorXd, Dynamic, Dynamic> E_;
	Array<MatrixXd, Dynamic, Dynamic> Theta_;

	Hmm() {}
	double backwardsWithCache(const Sequence& x, MatrixXd& betaBlock, MatrixXd& betaCache, VectorXd& betaLog2Scale);
	template<class Archive>
  void save(Archive & ar) const {
      const std::string type = "HMM";
      CEREALIZE(ar, type, Type);
      CEREALIZE(ar, nU_, nU); CEREALIZE(ar, nO_, nO); CEREALIZE(ar, dim_, dim);
      CEREALIZE(ar, T_, T); CEREALIZE(ar, E_, E); CEREALIZE(ar, pi_, pi);
  }
  template<class Archive>
  void load(Archive & ar) {
      std::string type = "HMM";
      CEREALIZE(ar, type, Type);
      CEREALIZE(ar, nU_, nU); CEREALIZE(ar, nO_, nO); CEREALIZE(ar, dim_, dim);
      CEREALIZE(ar, T_, T); CEREALIZE(ar, E_, E); CEREALIZE(ar, pi_, pi);
      init();
  }
};

/* ------------------------------------------------------------------------- *
 * Implementation                                                            *
 * ------------------------------------------------------------------------- */

void Hmm::setSize(int nStates, int nObservations, int nInputs, bool zeroParameters) {
	dim_ = nStates; nO_ = nObservations; nU_ = nInputs; nU1_ = (nU_ == 0 ? 1 : nU_);
	pi_.setConstant(dim_, zeroParameters ? 0 : double(1)/dim_);
	T_.setConstant( nU1_, MatrixXd::Constant(dim_, dim_, zeroParameters ? 0 : double(1)/dim_) );
	E_.setConstant( nO_, nU1_, VectorXd::Constant(dim_, zeroParameters ? 0 : double(1)/nO_) );
}

bool Hmm::normalize() {
	bool success = true;
	success &= (tom::normalize(pi_) > 0);
	for (int a = 0; a < nU1_; ++a) {
		success &= normalizeRows(T_(a));
		ArrayXd sum_Eo = E_.col(a).sum();
		if ((sum_Eo == 0).any()) { success = false; }
		else { for (int o = 0; o < nO_; ++o) { E_(o,a).array() /= sum_Eo; } }
	}
	return success;
}

void Hmm::randomize(double exponent, const Random& rnd) {
	Random& r = const_cast<Random&>(rnd);
	pi_ = r.random(dim_, 1).array().pow(exponent);
	for (int a = 0; a < nU1_; ++a) { T_(a) = r.random(dim_, dim_).array().pow(exponent); }
	for (int a = 0; a < nU1_; ++a) { for (int o = 0; o < nO_; ++o) { E_(o,a) = r.random(dim_,1).array().pow(exponent); }}
	normalize();
	init();
}

void Hmm::init() {
	nU1_ = (nU_ == 0 ? 1 : nU_);
	Theta_.resize(nO_, nU1_);
	for (int a = 0; a < nU1_; ++a) { for (int o = 0; o < nO_; ++o) { Theta_(o,a) = T_(a) * E_(o,a).asDiagonal(); }}
}

double Hmm::backwardsWithCache(const Sequence& x, MatrixXd& betaBlock, MatrixXd& betaCache, VectorXd& betaLog2Scale) {
	int N = x.length();
	int betaBlockSize = betaBlock.cols();
	VectorXd temp(dim_), beta(dim_);
	beta /*N-1*/ = VectorXd::Ones(dim_);
	betaLog2Scale(N-1) = 0;
	for (int t = N-1, tb = t / betaBlockSize, tbi = t % betaBlockSize; t >= betaBlockSize; --t) {
		temp.noalias() = Theta_(x.o(t), x.u(t)) * beta;
		beta /*t-1*/ = temp;
		betaLog2Scale(t-1) = betaLog2Scale(t) + std::log2( tom::normalize(beta) );
		if (tbi == 0) betaCache.col(tb-1) = beta;
		if (--tbi < 0) { tb--; tbi = betaBlockSize - 1; } // ensure t = tb * betaBlockSize + tbi
	}
	betaBlock.col(betaBlockSize-1) = betaCache.col(0);
	for (int tbc = betaBlockSize-1; tbc > 0; --tbc) {
		betaBlock.col(tbc-1).noalias() = Theta_(x.o(tbc), x.u(tbc)) * betaBlock.col(tbc);
		betaLog2Scale(tbc-1) = betaLog2Scale(tbc) + std::log2( tom::normalize( betaBlock.col(tbc-1) ) );
	}
	beta /*0*/ = betaBlock.col(0);

	if (nU_ != 0) { return betaLog2Scale(0) + std::log2( (beta.transpose() * Theta_(x.o(0), x.u(0)).transpose() * pi_).value() ); }
	else { return betaLog2Scale(0) + std::log2( (beta.transpose() * E_(x.o(0), x.u(0)).asDiagonal() * pi_).value() ); }
}

double Hmm::trainEM(const Sequence& trainSequence, const EMStopCondition& stopCondition) {
	int N = trainSequence.length(); if (N == 0) return 0;
	const Sequence& x = trainSequence;
	bool pomdp = (nU_ != 0);
	Hmm hmmOpt = Hmm(dim_, nO_, nU_, 0);
	double llOpt = 0;
	double log2Px;

	int betaBlockSize = ceil(sqrt(N));
	if (betaBlockSize < 3) { betaBlockSize = 1; }
	int betaCacheSize = (N-1) / betaBlockSize;
	MatrixXd betaCache(dim_, betaCacheSize);
	MatrixXd betaBlock(dim_, betaBlockSize);
	VectorXd betaLog2Scale(N);
	double alphaLog2Scale;

	VectorXd alpha(dim_);
	VectorXd beta(dim_);
	VectorXd gamma(dim_);
	VectorXd temp(dim_);
	MatrixXd xi(dim_, dim_);

	for (int it = 0;; ++it) { // The termination condition is checked in the loop
		log2Px = backwardsWithCache(x, betaBlock, betaCache, betaLog2Scale);
		llOpt = -log2Px / N;
		if (const_cast<EMStopCondition&>(stopCondition)(it, llOpt)) { return llOpt; }
		const_cast<EMStopCondition&>(stopCondition).previousLog2Likelihood_ = llOpt;

		hmmOpt.setSize(dim_, nO_, nU_, true);
		beta /*0*/ = betaBlock.col(0);
		// We first handle the starting case (t0 == -1 or t0 == 0)
		if (pomdp) {
			alpha /*-1*/ = pi_; // beta is still beta_0
			alphaLog2Scale = 0;
			beta /*-1*/ = Theta_(x.o(0), x.u(0)) * beta; tom::normalize(beta);
			hmmOpt.pi_ = alpha.cwiseProduct(beta); tom::normalize(hmmOpt.pi_); // = gamma_{-1}
			//alpha /*0*/ = Theta_(x.o(0), x.u(0)).transpose() * alpha; tom::normalize(alpha);
		} else { // HMM
			alpha /*0*/ = E_(x.o(0), x.u(0)).asDiagonal() * pi_;
			alphaLog2Scale = std::log2( tom::normalize(alpha) );
			gamma /*0*/ = alpha.cwiseProduct(beta); tom::normalize(gamma);
			hmmOpt.pi_ = gamma;
			hmmOpt.E_(x.o(0), x.u(0)) = gamma;
		}

		for (int t = ( pomdp ? 0 : 1 ), tb = t / betaBlockSize, tbi = t % betaBlockSize; t < N; ++t) {
			if (tbi == 0) { // recompute betaBlock
				if (tb == betaCacheSize) { // we are passed last cached value
					betaBlock.col( (N-1) % betaCacheSize ).setOnes();
					for (int tbc = (N-1) % betaCacheSize; tbc > 0; --tbc) {
						betaBlock.col(tbc-1).noalias() = Theta_(x.o(t+tbc), x.u(t+tbc)) * betaBlock.col(tbc);
						tom::normalize( betaBlock.col(tbc-1) );
					}
				} else {
					betaBlock.col(betaBlockSize-1) = betaCache.col(tb);
					for (int tbc = betaBlockSize-1; tbc > 0; --tbc) {
						betaBlock.col(tbc-1).noalias() = Theta_(x.o(t+tbc), x.u(t+tbc)) * betaBlock.col(tbc);
						tom::normalize( betaBlock.col(tbc-1) );
					}
				}
			}
			beta /*t*/ = betaBlock.col(tbi);
			xi.noalias() /*t-1*/ = std::exp2( alphaLog2Scale + betaLog2Scale(t) - log2Px ) * alpha.asDiagonal() * Theta_(x.o(t),x.u(t)) * beta.asDiagonal();
			// tom::normalize(xi);
			hmmOpt.T_(x.u(0)) += xi;
			temp.noalias() = Theta_(x.o(t), x.u(t)).transpose() * alpha;
			alpha /*t*/ = temp;
			alphaLog2Scale += std::log2( tom::normalize(alpha) );
			gamma /*t*/ = std::exp2( alphaLog2Scale + betaLog2Scale(t) - log2Px ) * alpha.cwiseProduct(beta);
			// tom::normalize(gamma);
			hmmOpt.E_(x.o(t), x.u(t)) += gamma;

			if (++tbi == betaBlockSize) { tb++; tbi = 0; } // ensure t = tb * betaBlockSize + tbi
		}

		hmmOpt.normalize();
		*this = hmmOpt;
		init();
	}
}

} // namespace tom

#endif // HMM_H
