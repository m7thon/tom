/**
 * @file   Oom.h
 * @author Michael Thon
 *
 * @brief  This file provides the basic functionality of observable operator models.
 */

#ifndef OOM_H
#define OOM_H

// SWIGCODE(%attribute(tom::Oom, int, nU, nU);)
SWIGCODE(%attribute(tom::Oom, int, maxSetback_, maxSetback, maxSetback);)
namespace tom {
SWIGCODE(%feature("python:slot", "tp_repr", functype="reprfunc") Oom::repr;)
/**
 * This class provides the basic functionality of observable operator models. The main aspects are
 * - extracting information from a given OOM (sequence generation, prediction, fingerprints, etc), and
 * - transformation functions (stabilization, equivalence transforms, dimension reduction, etc.)
 *  To estimate ("learn") an OOM from data, use the OomTrain class.
 */
class Oom {
	friend class cereal::access;
public:
/** @name Constructors and initialization*/ //@{
  /**
	 * Construct an uninitialized (!) Oom.
	 */
	Oom() {};
  /**
   * Construct a simple (random) OOM of dimension \a dim that models a stochastic process with \a nO number of possible observations and \a nU number of possible inputs.\ The initial state will be set to the stationary state in the case of an output-only \a Oom.
	 * @param nO the size of the output alphabet
	 * @param nU the size of the input alphabet (0 for an output-only \a Oom)
	 * @param dim the dimension of the OOM to construct
	 * @param exponent exponent of value distribution of tau operator matrix entries. When set to 1, the entries of matrices are sampled from a uniform distribution. Higher values lead to sample distributions that are increasingly skewed to the right, while a value of zero (default) will lead to an \a Oom that produces iud outputs.
   */
  Oom(int dim, int nO, int nU = 0, double exponent = 0, const Random& r = Random());
  /** Construct an \a Oom from the string representation given by \a oomString\. This must correspond to what the \a toString member function produces. */
  Oom(const std::string& oomString) { std::stringstream iss(oomString); *this << iss; }
	/** setup the internal structure for an OOM of the desired size without performing any initialization.\ Typically, the parameters \a sig, \a tau(o,u) and \a w0 will be assigned next, and then \a initialize() must be called.
	 */
	void setSize(int dim, int nO, int nU = 0);
	/** initialize the OOM.\ This assumes that all essential parameters (e.g., \a nU_, \a nO_, \a dim_, \a sig_, \a tau_, \a w0_ ) have been set.
	 */
  void init();

  void validate();
//@}

/** @name Accessors */ //@{
	/** return the size of the input alphabet */
	int nU() const { return nU_; }
	/** return the size of the output alphabet */
  int nO() const { return nO_; }
	/** return the model dimension */
	int dim() const { return dim_; }
  /** return a const reference to the evaluation functional vector \f$\sigma\f$ */
	const RowVectorXd& sig() const { return sig_; }
  /** return a reference to the evaluation functional vector \f$\sigma\f$ */
	RowVectorXd& sig() { return sig_; }
  /** set the evaluation functional vector \f$\sigma\f$ to \a s and reset() */
	void sig(const RowVectorXd& s) { sig_ = s;
		for (Symbol u = 0; u < (nU_ == 0 ? 1 : nU_); u++)
			for (Symbol o = 0; o < nO_; o++)
				sig_tau_(u).row(o) = sig_ * tau_(o, u);
		reset();
	}
	/** return a const reference to the observable operator corresponding to observation \a o and input \a u. */
	SWIGCODE(%feature("docstring", "return a const reference to the observable operator corresponding to observation \a o and input \a u.") tau;)
	const MatrixXd& tau(Symbol o, Symbol u = 0) const { return tau_(o, u); }
	/** return a reference to the observable operator corresponding to observation \a o and input \a u. */
	MatrixXd& tau(Symbol o, Symbol u = 0) { return tau_(o, u); }
	/** set the tau operator to the given matrix */
  void tau(Symbol o, Symbol u, const MatrixXd& tau_o_u) { tau(o,u) = tau_o_u; }
	/** return a const reference to the initial state */
	const VectorXd& w0() const { return w0_; }
	/** return a reference to the initial state */
	VectorXd& w0() { return w0_; }
	/** set the initial state to \a w and reset() */
	void w0(const VectorXd& w) { w0_ = w; reset();  }
	/** return a reference to the current state */
	const VectorXd& wt() { return wt_; }
	/** set the current state to \a w */
	void wt(const VectorXd& w, Sequence h = Sequence()) { history(h); didSetback_ = false; wt_ = w; if (nU_ == 0) condition(); }
	/** return a reference to the current prediction vector of the next output symbol probabilities, i.e., \f$P(\cdot|u_t, \omega_t)\f$. */
	const VectorXd& prediction() { return P_; }
	/** return the probability of the next output symbol \a o, i.e., \f$P(o|u_t, \omega_t)\f$. */
	double prediction(Symbol o) { return P_(o); }
	/** return the maximum number of steps to "replay" during a setBack() operation */
	int maxSetback() const { return maxSetback_; }
	/** set the maximum number of steps to "replay" during a setBack() operation to \a value */
	void maxSetback(int value) { maxSetback_ = value;
		out_buf_.resize(maxSetback_); in_buf_.clear(); in_buf_.resize(maxSetback_, 0); }
	/** set the history to the given \a Sequence, i.e., assume that the given \a Sequence has lead to the current state. */
	void history(const Sequence seq);
//@}

/** @name Basic OOM functionality */ //@{
  /** reset the OOM to its initial state and reset the error counters. */
  void reset() {
    historyLength_ = 0;
		didSetback_ = false;
    wt_ = w0_;
    if (nU_ == 0) condition();
    resetCounters();
  }
	/** reset the error counters */
	void resetCounters() { 	nSetback_ = 0; nFixPrediction_ = 0; nImpossible_ = 0; }
  /** update the Oom state according to the input-output pair (\a u,\a o), i.e., set \f$\omega_t = \tau_{u,o} \omega_t\f$, and \a normalize the state\. Then call \a condition() for OOMs, but not(!) for IO-OOMs\. This function performs setbacks when needed.*/
  void update(Symbol o, Symbol u = 0);
  /** update and normalize the prediction vector of the next output symbol probabilities according to the current state \f$\omega_t\f$ and input \a u, i.e., compute \f$P(\cdot|u_t, \omega_t)\f$\. This function calls setbacks when needed. */
  void condition(Symbol u = 0);
  /** generate a sample sequence of a given length according to the OOM and an input policy (by default an iud-process is used as input) starting from the current state \f$\omega_t\f$. */
  Sequence generate(unsigned long length, Random& r, const Policy& policy = Policy());
  /** generate a sample sequence of a given length according to the OOM and an input policy (by default an iud-process is used as input) starting from the current state \f$\omega_t\f$, but sample each observation from the probability vector to the power of \a alpha */
  Sequence generate(unsigned long length, Random& r, double alpha, const Policy& policy = Policy());
  /** return the function value (probability) for the input-output pair (\a u,\a o) given the current state, i.e., \f$P(o|u, \omega_t)\f$, and update the state. */
  double f(Symbol o, Symbol u = 0);
	/** return the function value (probability) for the given sequence \a seq given the current state, i.e., \f$f(seq | \omega_t)\f$, and update the state. */
	double f(const Sequence& seq);
	/** return the matrix \a F of function values \f$F^{I,J}=[f(\bar{x}_j\bar{x}_i)]_{i,j}\f$ for the given set \a chaSeqs of characteristic sequences \f$\bar{x}_j\f$ and the set \a indSeqs of indicative sequences \f$\bar{x}_j\f$. */
	MatrixXd f(const Sequences& chaSeqs, const Sequences& indSeqs);
	/** return the matrix \a Fz of function values \f$F_z^{I,J}=[f(\bar{x}_jz\bar{x}_i)]_{i,j}\f$ for the given set \a chaSeqs of characteristic sequences \f$\bar{x}_j\f$, the set \a indSeqs of indicative sequences \f$\bar{x}_j\f$ and the input-output pair \a z = (\a u,\a o). */
	MatrixXd f(const Sequences& chaSeqs, const Sequences& indSeqs, Symbol o, Symbol u = 0);
	/** return the log function value of the given sequence \a seq given the current state, i.e. \f$\log f(seq | \omega_t)\f$, and update the state.\ To deal gracefully with time-steps \a n in the sequence \a seq where \f$ P(o_n | u_n, \omega_{t+n}) \approx 0\f$, all such probabilities are treated as \a impossibleProbMargin_ in the computation of the function value.\ Every time this happens, the counter \a nImpossible_ is incremented.\ Note that this problem can be avoided by increasing \a minPredictionProb_. */
	double log_f(const Sequence& seq);

	/** return -log_2( f( \a seq ) ) / seq.length() after performing a \a reset() operation. */
	double ll(const Sequence& seq);
  /** calculate the average one-step squared prediction error along the given sample sequence \a seq according to a correct \a Oom model \a gen.\ Both this \a Oom and the given \a gen are evaluated from their current states, and their states are updated accordingly. */
	double averageOneStepPredictionError(Oom& gen, const Sequence& seq);
	/**
	 *
	 */
  MatrixXf* harvestStates(const Sequence& seq);
//@}

/** @name Transformation functions */ //@{
  /** return the stationary state (in the case of an input-output \a Oom according to a uniformly random iid input policy). */
	Eigen::VectorXd stationaryState() const;
  /** return the "reverse" OOM. */
  std::shared_ptr<Oom> reverse(int normalization = 1) const;
  /** transform the OOM to an equivalent OOM that has given \a sig and \a w0 as parameters for \a sig_ and \a w0_\. This will only yield an (equivalent) OOM if \a sig * \a w0 = 1. */
  void transform(const RowVectorXd& sig, const VectorXd& w0 = VectorXd::Zero(0));
	/** conjugate the OOM by the given matrices \a rho and \a rhoInv, i\.e\., set w0 = rho*w0, tau(o,u) = rho * tau(o,u) * rhoInv and sig = sig * rhoInv. **/
	void conjugate(const MatrixXd& rho, const MatrixXd& rhoInv);
	bool operator==(const Oom& other) const {
		bool equal = true;
		if (dim() != other.dim() or nO() != other.nO() or nU() != other.nU()) return false;
		if (sig() != other.sig() or w0() != other.w0()) return false;
		for (Symbol o = 0; o < nO(); ++o)
			for (Symbol u = 0; u < ((nU() == 0) ? 1 : nU()); ++u)
				if (tau(o,u) != other.tau(o,u)) return false;
		if (maxSetback() != other.maxSetback() or
				minPrediction_ != other.minPrediction_ or
				maxPredictionError_ != other.maxPredictionError_) return false;
		return true;
	}
//@}

/** @name IO-functions */ //@{
  /** return a \a std::string representation that can be used for saving to file etc. This function should just call the output stream operator. */
	std::string toString() const { std::stringstream oss; *this >> oss; return oss.str(); }
	INSERT_JSON_IO_FUNCTIONS()
	/** return a representation to display in interactive python. */
	std::string repr() const;
  /** read from the given input stream and initialize\. The format must correspond to what the output function produces. */
	std::istream & operator<<(std::istream &istream);
  /** write to the given output stream. */
	std::ostream& operator>>(std::ostream &os) const;
//@}

/** @name Internalals\. Use only if you know what you are doing! */ //@{
  /** attempt to fix the prediction vector of the next output symbol probabilities \f$P(\cdot|u_t, \omega_t)\f$ such that all probabilities are at least \a minPredictionProb_ and the probabilities sum to one.\ Return a measure of the required change to the prediction vector: 1.5 * nO() * squared norm of the difference. */
	double normalizePrediction();
	/** attempt to perform a state setback operation for at most \a maxSetback_ time-steps.\ Note that calling this method repeatedly will attempt a setback for a shorter history each time.\ Return \c true if a setback could be performed. */
	bool setBack();
//@}

	bool valid_ =              false; /**< if this OOM is guaranteed to be valid, we may be able to skip some stabilization heuristics */
  double minPrediction_ = 0 /*1e-5*/; /**< the smallest probability to assign to any next-symbol prediction when normalizing predictions */
  double maxPredictionError_ =   0.3; /**< the largest tolerated error of the prediction vector before a state setback is performed */
	double fixPredictionMargin_ = 0.01; /**< the largest tolerated error of the prediction vector before the normalization is considered a fixing event */
  int nSetback_ =                  0; /**< the number of times that the OOM state needed to be fixed by a setback operation (introduces a considerable error) */
  int nFixPrediction_ =            0; /**< the number of times that the predicted probabilities were adjusted by more than \a fixPredictionMargin_. */
	int nImpossible_ =               0; /**< the number of times that a symbol was encountered that has probability smaller than \a impossibleProbMargin_ according to the OOM */
	double epsilonZero_ =        1e-12; /**< the smallest number reasonably considered non-zero (e.g., for division purposes) */
	int outputPrecision_ = std::numeric_limits< double >::digits10; /**< the precision to use when outputting the oom */

private:
	int nU_;           /**< the size of the input alphabet */
  int nO_;           /**< the size of the output alphabet */
  int dim_;          /**< the model dimension */
	RowVectorXd sig_;  /**< the probability functional vector \f$\sigma\f$ */
	Array<MatrixXd, Dynamic, Dynamic> tau_; /**< the array of observable operators (of size \a nO x \a nU). The operator corresponding to the input-observation pair (u,o) is given by tau(o,u). (Read: "the observable operator for the observation o given the input u"). This is convenient when dealing with standard OOMs (\a nU == 0), since one can address the operators as tau(o) in this case, as the action u defaults to 0). */
	VectorXd w0_;      /**< the initial state */
	VectorXd wt_;      /**< the current state (at time t) */
  unsigned int maxSetback_ = 0 /*3*/; /**< the maximum number of steps to "replay" during a setBack operation */

  bool didSetback_ = false;         /**< \c true if a setback operation has been performed in the current time step */
	Array<MatrixXd, 1, Dynamic> sig_tau_;
  std::deque<Symbol> in_buf_, out_buf_;
	VectorXd P_;       /**< the conditional next-symbol probabilities given the current input: \f$P(o)=P(o|u_t, w_t)=\sigma*tau(o,u)*w_t\f$. */
  unsigned long historyLength_ = 0;
	VectorXd temp_dim_, temp_nO_;

  /** a helper function to read optional parameters from stream */
  bool readOptionalPropertyFromStream(std::istream &istream);

  /** a helper function to write optional parameters to stream */
  void writeOptionalPropertiesToStream(std::ostream &ostream) const;

#ifndef SWIG
	template<class Archive>
  void save(Archive & ar) const {
		const std::string type = "OOM";
		ar(cereal::make_nvp("Type", type));
		MVAR(ar,nU); MVAR(ar,nO); MVAR(ar,dim);
		MVAR(ar,sig); MVAR(ar,tau); MVAR(ar,w0);
		MVAR(ar,minPrediction); MVAR(ar,maxPredictionError); MVAR(ar,maxSetback);
  }

  template<class Archive>
  void load(Archive & ar) {
		std::string type;
		ar(cereal::make_nvp("Type", type));
		MVAR(ar,nU); MVAR(ar,nO); MVAR(ar,dim);
		MVAR(ar,sig); MVAR(ar,tau); MVAR(ar,w0);
		OMVAR(ar,minPrediction); OMVAR(ar,maxPredictionError); OMVAR(ar,maxSetback);
		init();
		validate();
	}
#endif // SWIG
};

#ifdef SWIG
%pythoncode %{
	OomNoStabilizationParams = { 'minPrediction_': 0, 'maxPredictionError_': 0.3, 'maxSetback_': 0, 'epsilonZero_': 1e-12 }
	OomDefaultStabilizationParams = { 'minPrediction_': 0.0003, 'maxPredictionError_': 0.02, 'maxSetback_': 3, 'epsilonZero_': 1e-8 }
%}
#endif

#ifndef SWIG
/**
 * write the parameters of the OOM \a oom to the given output stream.
 */
std::ostream& operator<<(std::ostream& os, const Oom& oom) { return oom >> os; }
/**
 * read the parameters of the OOM \a oom from the given input stream and initialize it. The format must correspond to what the output functions produce.
 */
std::istream& operator>>(std::istream& istream, Oom& oom) { return oom << istream; }
#endif // SWIG

} // namespace tom

#endif // OOM_H
