#include "tom.h"

namespace tom {

// CONSTRUCTORS AND INITIALIZATION
Oom::Oom(int dim, int nO, int nU, double exponent, const Random& r) : dim_(dim), nO_(nO), nU_(nU) {
	Random& rnd = const_cast<Random&>(r);
	setSize(dim, nO, nU);
	if (exponent == 0) return;
	for (Symbol u = 0; u < (nU_ == 0 ? 1 : nU_); u++) {
		for (Symbol o = 0; o < nO_; o++) {
			for (int i = 0; i < dim_; i++) {
				for (int j = 0; j < dim_; j++) {
					tau_(o, u)(i,j) = std::pow(rnd.random(), exponent);
				}
			}
		}
		MatrixXd tau = MatrixXd::Zero(dim_,dim_);
		for (Symbol o = 0; o < nO_; o++) {
			tau += tau_(o, u);
		}
		RowVectorXd sig_tau = tau.colwise().sum();
		for (Symbol o = 0; o < nO_; o++) {
			for (int j = 0; j < dim_; j++) {
				tau_(o,u).col(j) /= sig_tau(j);
			}
		}
	}
	w0_ = stationaryState();
  init();
}

void Oom::setSize(int dim, int nO, int nU) {
	nU_ = nU; nO_ = nO, dim_ = dim;
	tau_.resize(nO_, (nU_ == 0 ? 1 : nU_));
	w0_ = VectorXd::Ones(dim_) / dim_;
	sig_ = RowVectorXd::Ones(dim_);
	for (Symbol u = 0; u < (nU_ == 0 ? 1 : nU_); u++) {
		for (Symbol o = 0; o < nO_; o++) {
			tau_(o, u) = MatrixXd::Identity(dim_, dim_) / nO_;
		}
	}
	temp_dim_.resize(dim_);
	temp_nO_.resize(nO_);
  init();
}

void Oom::init() {
	maxSetback(maxSetback_);
	sig_tau_.resize(1, (nU_ == 0 ? 1 : nU_));
	P_.resize(nO_, 1);
	for (Symbol u = 0; u < (nU_ == 0 ? 1 : nU_); u++) {
		sig_tau_(u).resize(nO_, dim_);
		for (Symbol o = 0; o < nO_; o++)
			sig_tau_(u).row(o) = sig_ * tau_(o, u);
	}
  reset();
}

void Oom::validate() {
  valid_ = true;
  if (((sig_ * w0_).array() - 1.0).pow(2)(0,0) > epsilonZero_)
    valid_ = false;
  for (Symbol u = 0; u < (nU_ == 0 ? 1 : nU_); u++)
    if ((sig_tau_(u).colwise().sum() - sig_).norm() > epsilonZero_)
      valid_ = false;
  // Next check if we are dealing with a POMDP
  if (w0_.minCoeff() < 0)
    valid_ = false;
  if (sig_.minCoeff() < 0)
    valid_ = false;
  const Symbol u_max = (nU_ == 0 ? 1 : nU_);
  for (Symbol u = 0; u < u_max; ++u)
    for (Symbol o = 0; o < nO_; ++o)
      if (tau_(o, u).minCoeff() < 0)
        valid_ = false;
}

// IO-FUNCTIONS
void Oom::writeOptionalPropertiesToStream(std::ostream &os) const {
  os << "minPrediction: " << minPrediction_ << std::endl;
  os << "maxPredictionError: " << maxPredictionError_ << std::endl;
  os << "maxSetback: " << maxSetback_ << std::endl;
  return;
}
std::ostream& Oom::operator>>(std::ostream &os) const {
	os << std::setprecision(outputPrecision_);
	os << "OOM" << " valid: " << (valid_ ? "True" : "False") << " nU: " << nU_ << " nO: " << nO_ << " dim: " << dim_
     << std::endl;
  os << "sig:" << std::endl << sig_ << std::endl;
	for (Symbol u = 0; u < (nU_ == 0 ? 1 : nU_); u++)
		for (Symbol o = 0; o < nO_; o++) {
			os << "tau(" << o;
			if (nU_ == 0) os << "):" << std::endl << tau_(o, u) << std::endl;
			else os << "," << u << "):" << std::endl << tau_(o, u) << std::endl;
		}
	os << "w0:" << std::endl << w0_ << std::endl;
  writeOptionalPropertiesToStream(os);
  os << "OOM_END" << std::endl;
	return os;
}
bool Oom::readOptionalPropertyFromStream(std::istream &istream) {
  std::string description; istream >> description;
  if (description == "OOM_END") { return false; }
  else if (description == "minPrediction:") { istream >> minPrediction_; }
  else if (description == "maxPredictionError:") { istream >> maxPredictionError_; }
  else if (description == "maxSetback:") { istream >> maxSetback_; }
  else {
    std::cerr << "unrecognized property: " << description << std::endl;
    istream >> description;
  }
  return true;
}
std::istream& Oom::operator<<(std::istream &istream) {
	std::string description;
  istream >> description;
  istream >> description >> description; valid_ = (description == "True");
	istream >> description >> nU_ >> description >> nO_ >> description >> dim_;
	w0_.resize(dim_);
	sig_.resize(dim_);
	tau_.resize(nO_, (nU_ == 0 ? 1 : nU_));
	istream >> description;
	for (int i = 0; i < dim_; ++i)
		istream >> sig_(i);
	for (Symbol u = 0; u < (nU_ == 0 ? 1 : nU_); u++)
		for (Symbol o = 0; o < nO_; o++)
		{
			istream >> description;
			tau_(o,u).resize(dim_,dim_);
			for (int i = 0; i < dim_; ++i)
				for (int j = 0; j < dim_; ++j)
					istream >> tau_(o,u)(i,j);
		}
	istream >> description;
	for (int i = 0; i < dim_; ++i)
		istream >> w0_(i);
  init();
	validate();
  while (readOptionalPropertyFromStream(istream));
	return istream;
}

std::string Oom::repr() const {
	std::stringstream os;
	os << "OOM" << " nU: " << nU_ << " nO: " << nO_ << " dim: " << dim_; return os.str();
}

// ACCESSORS
void Oom::history(const Sequence seq) {
	historyLength_ = seq.length();
	for (int i = ( historyLength_ < maxSetback_ ? historyLength_ : maxSetback_ ); i > 0; i--) {
		out_buf_[maxSetback_ - i] = seq.o(historyLength_ - i);
		if (nU_ != 0) { in_buf_[maxSetback_ - i] = seq.u(historyLength_ - i); }
	}
}

// CORE FUNCTIONS
void Oom::update(Symbol o, Symbol u) {
  // update state
	temp_dim_.noalias() = tau_(o,u) * wt_; wt_ = temp_dim_;
	out_buf_.push_back(o); out_buf_.pop_front(); ++historyLength_;
	if (nU_ != 0) { in_buf_.push_back(u); in_buf_.pop_front(); }
  didSetback_ = false;
  // normalize state
  double sig_wt = sig_ * wt_;
	if (sig_wt < epsilonZero_) setBack();
  else wt_ /= sig_wt;
  // if not IO, condition, i.e., get and fix probability vector:
  if (nU_ == 0) condition();
}
double Oom::normalizePrediction() {
	temp_nO_ = P_;
	P_.array() -= minPrediction_;
	P_ = P_.cwiseMax(0);
	double P_sum = P_.sum();
	if (P_sum < epsilonZero_) P_.setConstant(1./nO_);
	else P_ = P_.array() * (1 - nO_ * minPrediction_)/P_sum + minPrediction_;
	return sqrt((P_ - temp_nO_).squaredNorm() / nO_);
	// This used to be (P_ - oldP_).squaredNorm() * 1.5 * nO_;
}
void Oom::condition(Symbol u) {
  P_.noalias() = sig_tau_(u) * wt_;
	double fP = 0;
	while ( ( ( fP = normalizePrediction() ) > maxPredictionError_ ) and setBack() )
    P_.noalias() = sig_tau_(u) * wt_;
	if (fP > fixPredictionMargin_) ++nFixPrediction_;
}
bool Oom::setBack() {
  if (historyLength_ == 0) { wt_ = w0_; return false; } // nothing to setback!
	if (!didSetback_) { nSetback_++; didSetback_ = true; }
  double sig_wt = 0;
  historyLength_ = std::min((unsigned long)(maxSetback_ + 1), historyLength_); // --historylength next...
  while ((sig_wt < epsilonZero_) and (--historyLength_ > 0)) {
    // we decrease history by one, otherwise we just get the same state back
    wt_ = w0_;
    for ( unsigned long i = historyLength_; i > 0; i-- ) {
      // replay the history
      temp_dim_ = tau_(out_buf_[out_buf_.size() - i], in_buf_[in_buf_.size() - i]) * wt_; wt_ = temp_dim_;
      sig_wt = sig_ * wt_;
      if (sig_wt < epsilonZero_) break; // still bad states: try with shorter history...
      wt_ /= sig_wt;
    }
  }
  if (historyLength_ == 0) wt_ = w0_;
  return true;
}

Sequence Oom::generate(unsigned long length, Random& r, const Policy& policy) {
  Sequence seq(length, nO_, nU_);
	Symbol u,o;
	for (unsigned long l = 0; l < length; l++) {
		LOOP_PROGRESS("Generating sequence", l, length)
		if (nU_ != 0) { // IO-OOM
			u = (policy.nU_ == 0 ? r.integer(nU()) : policy.u(wt(), r));
			condition(u);
			o = r.sample(P_);
			seq.u(l, u); seq.o(l, o);
			update(o, u);
		}
		else { // OOM
			o = r.sample(P_);
			seq.o(l, o);
			update(o);
		}
	}
	LOOP_DONE("Generating sequence")
	return seq;
}
Sequence Oom::generate(unsigned long length, Random& r, double alpha, const Policy& policy) {
  Sequence seq(length, nO_, nU_);
	Symbol u,o;
	for (unsigned long l = 0; l < length; l++) {
		LOOP_PROGRESS("Generating sequence", l, length)
		if (nU_ != 0) { // IO-OOM
			u = (policy.nU_ == 0 ? r.integer(nU()) : policy.u(wt(), r));
			condition(u);
			o = r.sample(P_.array().pow(alpha).matrix());
			seq.u(l, u); seq.o(l, o);
			update(o, u);
		}
		else { // OOM
			o = r.sample(P_.array().pow(alpha).matrix());
			seq.o(l, o);
			update(o);
		}
	}
	LOOP_DONE("Generating sequence")
	return seq;
}

double Oom::f(Symbol o, Symbol u) {
	if (nU_ != 0) condition(u);
  double val = prediction(o);
	update(o, u);
	return val;
}
double Oom::f(const Sequence& seq) {
	double val = 1;
	for (unsigned long pos = 0; pos < seq.length(); ++pos)
		val *= f( seq.o(pos), ((nU_ != 0) ? seq.u(pos) : 0) );
	return val;
}
MatrixXd Oom::f(const Sequences& chaSeqs, const Sequences& indSeqs) {
	unsigned int Frows = chaSeqs.size(), Fcols = indSeqs.size();
	MatrixXd F = MatrixXd::Zero(Frows, Fcols);
	for (unsigned int j = 0; j < Fcols; ++j) {
		reset();
		double f_temp = f(indSeqs[j]);
		VectorXd wt_temp = wt_;
		for (unsigned int i = 0; i < Frows; ++i) {
			wt(wt_temp, indSeqs[j]);
			F.coeffRef(i,j) = f_temp * f(chaSeqs[i]);
		}
	}
	return F;
}
MatrixXd Oom::f(const Sequences& chaSeqs, const Sequences& indSeqs, Symbol o, Symbol u) {
	unsigned int Frows = chaSeqs.size(), Fcols = indSeqs.size();
	MatrixXd F = MatrixXd::Zero(Frows, Fcols);
	for (unsigned int j = 0; j < Fcols; ++j) {
		reset();
		double f_temp = f(indSeqs[j]);
		VectorXd wt_temp = wt_;
		for (unsigned int i = 0; i < Frows; ++i) {
			wt(wt_temp, indSeqs[j]);
			F.coeffRef(i,j) = f_temp * f(o,u) * f(chaSeqs[i]);
		}
	}
	return F;
}

double Oom::log_f(const Sequence& seq) {
	double LL = 0;
	double val;
	for (unsigned long pos = 0; pos < seq.length(); ++pos) {
		val = f( seq.o(pos), seq.u(pos) );
		if (val < epsilonZero_) {
      LL += log(epsilonZero_);
			nImpossible_++;
		}	else {
			LL += log(val);
		}
	}
	return LL;
}

double Oom::ll(const Sequence& seq) {
	reset();
	static const double log2_e = std::log2(std::exp(1.0));
	return -log2_e * log_f(seq) / seq.length();

}

double Oom::averageOneStepPredictionError(Oom& gen, const Sequence& seq) {
	double sqr_err = 0;
	for (unsigned long pos = 0; pos < seq.length(); pos++) {
		Symbol u = seq.u(pos); Symbol o = seq.o(pos);
		condition(u); gen.condition(u);
		sqr_err += (prediction() - gen.prediction()).squaredNorm();
		update(o, u); gen.update(o, u);
	}
	return sqr_err / ( seq.length() * nO_ );
}

MatrixXf* Oom::harvestStates(const Sequence& seq) {
  double LLRevTrain = 0;
  long l = seq.length();

  MatrixXf* W = new MatrixXf(dim_, l+1);

  W->col(0) = wt().cast<float>();
  for (long i = 0; i < l; ++i) {
    double val = f(seq.o(i), ((nU_ != 0) ? seq.u(i) : 0));
    W->col(i+1) = wt().cast<float>();
    if (val < epsilonZero_) {
      LLRevTrain += log(epsilonZero_);
      nImpossible_++;
    }	else {
      LLRevTrain += log(val);
    }
  }
  return W;
}

// TRANSFORMATIONS
VectorXd Oom::stationaryState() const {
	MatrixXd t = MatrixXd::Zero(dim_,dim_);
  for (Symbol u = 0; u < (nU() == 0 ? 1 : nU()); u++) {
    for (Symbol o = 0; o < nO_; o++) { t += tau(o, u); }
	}
	t /= (nU() == 0 ? 1 : nU());
	VectorXcd eigs; MatrixXcd eigv;
	eigensolve(eigs, eigv, t);
	// find eigenvalue closest to 1
	double distFromOne = std::abs(eigs(0) - 1.);
	int oneIndex = 0;
	for (int i = 1; i < eigs.size(); ++i) {
		if (std::abs(eigs(i) - 1.) < distFromOne) {
			distFromOne = std::abs(1. - eigs(i));
			oneIndex = i;
		}
	}
	VectorXd wStationary = eigv.col(oneIndex).real();
	wStationary /= (double)(sig() * wStationary);
  return wStationary;
}

std::shared_ptr<Oom> Oom::reverse(int normalization) const {
  std::shared_ptr<Oom> room(new Oom(*this)); // copy of oom;
  MatrixXd rho = MatrixXd::Identity(dim(), dim());
  MatrixXd rhoInv = MatrixXd::Identity(dim(), dim());
  if (normalization == 1) {
    rho = w0().asDiagonal().inverse(); // rho * w0 = (1,...,1)
    rho /= dim(); // it should also work without this
    rhoInv = dim() * w0().asDiagonal();
  }
  for (Symbol a = 0; a < (nU_ == 0 ? 1 : nU_); ++a)
    for (Symbol o = 0; o < nO(); ++o)
      room->tau(o, a) = (rho * tau(o, a) * rhoInv).transpose();
  room->w0((sig() * rhoInv).transpose());
  room->sig((rho * w0()).transpose());
  room->reset();
  return room;
}

void Oom::conjugate(const MatrixXd& rho, const MatrixXd& rhoInv) {
	w0() = rho * w0();
	for (Symbol u = 0; u < (nU_ == 0 ? 1 : nU_); u++)
		for (Symbol o = 0; o < nO_; o++)
			tau(o, u) = rho * tau(o,u) * rhoInv;
	sig() = sig() * rhoInv;
	dim_ = rho.rows();
	init();
}

void Oom::transform(const RowVectorXd& sig, const VectorXd& w0) {
	VectorXd e1 = VectorXd::Zero(dim()); e1(0) = 1;
	VectorXd v = sig_.transpose() / sig_.norm() - e1;
	MatrixXd H = MatrixXd::Identity(dim(), dim());
	if (not v.isZero()) {
		v /= v.norm();
		H -= 2 * v * v.transpose();
	}
	H.col(0) = w0_;
	v = sig.transpose() / sig.norm() - e1;
	MatrixXd H2 = MatrixXd::Identity(dim(), dim());
	if (not v.isZero()) {
		v /= v.norm();
		H2 -= 2 * v * v.transpose();
	}
	H2.col(0) = w0;
	conjugate(H2 * H.inverse(), H * H2.inverse());
}

} // namespace tom
