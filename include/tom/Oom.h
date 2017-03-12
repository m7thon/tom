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
    /* <editor-fold desc="Constructors and initialization"> */ /** @name Constructors and initialization */ //@{
    /**
     * Construct an uninitialized (!) Oom. */
    Oom() { stabilization(-1, -1, -1, -1, "none"); };

    /** Construct a simple (random) OOM of dimension \c dimension that models a stochastic process with \c nOutputSymbols number of possible observations and \c nInputSymbols number of possible inputs. The initial state will be set to the stationary state, assuming iid inputs in the case of an input-output OOM.
     * @param dimension the dimension of the OOM
     * @param nOutputSymbols the size of the output alphabet
     * @param nInputSymbols the size of the input alphabet, or 0 (default) for an output-only \c Oom
     * @param randomExponent exponent of value distribution of tau operator matrix entries. When set to 1, the entries of matrices are sampled from a uniform distribution. Higher values lead to sample distributions that are increasingly skewed to the right, while a value of zero (default) will lead to an \c Oom that produces iid outputs.
     * @param zero_threshold set all parameters less than `zero_threshold` to zero and renormalize.
     * @param randomSource the `Random` object to use as the source of randomness
     */
    Oom(int dimension, int nOutputSymbols, int nInputSymbols = 0, double randomExponent = 0, double zero_threshold = 0, const Random &randomSource = Random()) :
            dim_(dimension), nO_(nOutputSymbols), nU_(nInputSymbols) {
        setSize(dimension, nOutputSymbols, nInputSymbols);
        if (randomExponent == 0) return;
        for (Symbol u = 0; u < (nU_ == 0 ? 1 : nU_); u++) {
            for (Symbol o = 0; o < nO_; o++)
                for (int i = 0; i < dim_; i++)
                    for (int j = 0; j < dim_; j++)
                        tau_(o, u)(i, j) = std::pow(randomSource.random(), randomExponent);
            MatrixXd tau = MatrixXd::Zero(dim_, dim_);
            for (Symbol o = 0; o < nO_; o++)
                tau += tau_(o, u);
            RowVectorXd sig_tau = tau.colwise().sum();
            for (Symbol o = 0; o < nO_; o++)
                for (int j = 0; j < dim_; j++)
                    tau_(o, u).col(j) /= sig_tau(j);
            if (zero_threshold > 0) {
                for (Symbol o = 0; o < nO_; o++)
                    for (int i = 0; i < dim_; i++)
                        for (int j = 0; j < dim_; j++) {
                            if (tau_(o, u)(i, j) < zero_threshold) tau_(o, u)(i, j) = 0;
                        }
                tau = MatrixXd::Zero(dim_, dim_);
                for (Symbol o = 0; o < nO_; o++)
                    tau += tau_(o, u);
                sig_tau = tau.colwise().sum();
                for (Symbol o = 0; o < nO_; o++)
                    for (int j = 0; j < dim_; j++)
                        tau_(o, u).col(j) /= sig_tau(j);
            }
        }
        w0_ = stationaryState();
        stabilization(-1, -1, -1, -1, "none");
        initialize();
    }

    /**
     * Construct an `Oom` corresponding to the given string `json_representation`. The format must correspond to what `toJSON()` produces. */
    Oom(const std::string &json_representation) { fromJSON(json_representation); }

    /**
     * Construct an `Oom` equivalent to the `Hmm` given by `hmm`. */
    Oom(const Hmm& hmm) {
        setSize(hmm.nStates(), hmm.nObservations(), hmm.nInputs());
        for (Symbol u = 0; u < (nU_ == 0 ? 1 : nU_); u++) {
            for (Symbol o = 0; o < nO_; o++) {
                if (isIO()) tau_(o,u) = hmm.E(o,u).asDiagonal() * hmm.T(u).transpose();
                else          tau_(o) = hmm.T().transpose() * hmm.E(o).asDiagonal();
            }
        }
        w0_ = hmm.pi();
        stabilization(-1, -1, -1, -1, "none");
        initialize();
    }

    /** Set the internal structure for an OOM of the desired size without performing any initialization.\ Typically, the parameters \c sig, \c tau(o,u) and \c w0 will be assigned next, and then \c initialize() must be called.
     * @param dimension the dimension of the OOM
     * @param nOutputSymbols the size of the output alphabet
     * @param nInputSymbols the size of the input alphabet, or 0 (default) for an output-only \c Oom
     */
    void setSize(int dimension, int nOutputSymbols, int nInputSymbols = 0) {
        dim_ = dimension;
        nO_ = nOutputSymbols;
        nU_ = nInputSymbols;
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
        initialize();
    }

    /** Initialize the OOM.\ This assumes that all essential parameters (i.e., `dimension`, `nOutputSymbols`, `nInputSymbols`, `sig`, `tau(o,u)` and `w0`) have been set.
     */
    void initialize() {
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
    //@} </editor-fold>

    /* <editor-fold desc="Main Oom parameters"> */ /** @name Main Oom parameters */ //@{
    /**
     * Return the model dimension */
    int dimension() const { return dim_; }

    /**
     * Return the size of the input alphabet\. Use `setSize()` to modify. */
    int nInputSymbols() const { return nU_; }

    /**
     * Return the size of the output alphabet\. Use `setSize()` to modify */
    int nOutputSymbols() const { return nO_; }

    /**
     * Return `true` if this is an input-output sequence, i.e., if the input alphabet size `nInputSymbols()` is non-zero. */
    bool isIO() const { return ( nU_ != 0 ); }

    /** Set this to be an IO-OOM if `io` is true, else to be an OOM. This only makes sense for OOMs/IO-OOMs with at most one input symbol. */
    void setIO(bool io = true) throw(std::runtime_error) {
        if (nU_ > 1) {
            if (io) return;
            throw std::runtime_error("Too many input symbols to convert to normal OOM!");
        }
        nU_ = (io ? 1 : 0);
    }

    /**
     * Return the evaluation functional row vector \f$\sigma\f$ */
    const RowVectorXd &sig() const { return sig_; }

    /**
     * Set the evaluation functional vector \f$\sigma\f$ to the given row vector \c new_value. Note that you *must* call `initialize()` after re-setting the `Oom` parameters `sig` or `tau(o,u)`. */
    void sig(const RowVectorXd &new_value) { sig_ = new_value; }

    /**
     * Return the observable operator corresponding to observation \c o and input \c u. The parameter \c u defaults to 0 for the case of no inputs. */
    const MatrixXd &tau(Symbol o, Symbol u = 0) const { return tau_(o, u); }

    /**
     * Return the observable operator corresponding to the symbol `z` given as a `Sequence` of length one. */
    const MatrixXd &tau(const Sequence& z) const throw(std::invalid_argument) {
        if (z.length() != 1 or nO_ != z.nOutputSymbols() or nU_ != z.nInputSymbols()) throw std::invalid_argument("invalid symbol z for this Oom");
        return tau_(z.o(0), z.u(0));
    }

    /**
     * Set the observable operator corresponding to the symbol `z` given as a `Sequence` of length one to the given matrix `new_value`. Note that you *must* call `initialize()` after re-setting the `Oom` parameters `sig` or `tau(o,u)`. */
    void tau(const Sequence& z, const MatrixXd &new_value) throw(std::invalid_argument) { const_cast<MatrixXd&>(tau(z)) = new_value; }

    /**
     * Set the observable operator corresponding to observation \c o and input \c u to the given matrix `new_value`. Note that you *must* call `initialize()` after re-setting the `Oom` parameters `sig` or `tau(o,u)`. */
    void tau(Symbol o, Symbol u, const MatrixXd &new_value) { tau_(o, u) = new_value; }

    /**
     * Set the observable operator corresponding to the observation \c o to the given matrix `new_value` */
    void tau(Symbol o, const MatrixXd &new_value) { tau(o, 0, new_value); }

    /**
     * Return the initial state vector \f$\omega_0\f$ */
    const VectorXd &w0() const { return w0_; }

    /**
     * Set the initial state vector \f$\omega_0\f$ to the given vector \c new_value and perform a \c reset(). */
    void w0(const VectorXd &new_value) {
        w0_ = new_value;
        reset();
    }
    //@} </editor-fold>

    /* <editor-fold desc="Stabilization"> */ /** @name Stabilization */ //@{
    SWIGCODE(%feature ("kwargs") stabilization;)
    SWIGCODE(%apply double *OUTPUT { double *minPredictionOut, double *normalizationToleranceOut, double *impossibilityThresholdOut };)
    SWIGCODE(%apply int *OUTPUT { int *maxSetbackOut };)
    /**
    Set (optional) and then return the stabilization parameters for the `Oom`. This is done as follows:
    - First, if a `preset` ("none" / "default") is specified, all stabilization parameters are set accordingly.
    - Next, any non-default argument causes the corresponding stabilization parameter to be set to the given value, while any argument left at its default value (-1) has no effect.
    - Finally, the current stabilization parameters are returned
    \if PY
    in a tuple in the same order as they appear as function arguments. This allows writing python code such as:
    \verbatim
    old_params = oom.stabilization()
    ...
    oom.stabilization(*old_params)
    \endverbatim
    \else
    in the corresponding output arguments `...Out` (if not NULL).
    \endif

    @param minPrediction The minimum probability for any observation symbol\. The `prediction()` is normalized at each time-step such that every output symbol has at least this probability.
    @param normalizationTolerance The maximum allowed tolerance between the unnormalized and normalized prediction vector before invoking a `setBack()`\. The tolerance is computed as 1.5 * `nOutputSymbols()` * squared norm of the difference\. See also `normalizePrediction()`.
    @param maxSetback The maximum number of steps to "replay" during a `setBack()` operation.
    @param impossibilityThreshold The smallest probability value considered as non-zero\. If an observation symbol is encountered that this model predicts to have a probability lower than `impossibilityThreshold()`, this is considered as an impossible event and a `setBack()` must be performed, as the state `wt()` can no longer be normalized.
    @param preset A preset for the stabilization parameters to apply before setting the other parameters.\ This can be either "none" to disable stabilization, or "default" to use the default stabilization settings.
    */
    C1(void) PY1(tuple)
    stabilization(C5(double *minPredictionOut, double *normalizationToleranceOut, int *maxSetbackOut, double *impossibilityThresholdOut,)
                  double minPrediction = -1, double normalizationTolerance = -1, int maxSetback = -1, double impossibilityThreshold = -1,
                  std::string preset = "") throw(std::invalid_argument) {
        if (preset == "default") {
            minPrediction_ = 0.0002;
            normalizationTolerance_ = 0.03;
            impossibilityThreshold_ = 1e-8;
            this->maxSetback(5);
        } else if (preset == "none") {
            minPrediction_ = 0;
            normalizationTolerance_ = std::numeric_limits<double>::infinity();
            impossibilityThreshold_ = 1e-15;
            this->maxSetback(0);
        } else if (preset != "") {
            throw std::invalid_argument("unrecognized preset");
        }
        if (minPrediction != -1) this->minPrediction(minPrediction);
        if (normalizationTolerance != -1) this->normalizationTolerance(normalizationTolerance);
        if (impossibilityThreshold != -1) this->impossibilityThreshold(impossibilityThreshold);
        if (maxSetback >= 0) this->maxSetback((unsigned int)maxSetback);
        if (minPredictionOut) *minPredictionOut = this->minPrediction_;
        if (normalizationToleranceOut) *normalizationToleranceOut = this->normalizationTolerance_;
        if(maxSetbackOut) *maxSetbackOut = this->maxSetback_;
        if(impossibilityThresholdOut) *impossibilityThresholdOut = this->impossibilityThreshold_;
    }
    SWIGCODE(%clear double *minPredictionOut, double *normalizationToleranceOut, double *impossibilityThresholdOut;)
    SWIGCODE(%clear int *maxSetbackOut;)

#ifndef SWIG
#ifndef PY
    /**
     * Set the stabilization parameters for the `Oom`. This is a convenience function that simply calls `stabilization(NULL, NULL, NULL, NULL, minPrediction, normalizationTolerance, maxSetback, impossibilityThreshold, preset);` */
    void stabilization(double minPrediction = -1, double normalizationTolerance = -1, int maxSetback = -1,
                       double impossibilityThreshold = -1, std::string preset = "") throw(std::invalid_argument) {
        stabilization(NULL, NULL, NULL, NULL, minPrediction, normalizationTolerance, maxSetback, impossibilityThreshold, preset);
    }
#endif
#endif

    /**
     * Return the minimum probability for any observation symbol. The `prediction()` is normalized at each time-step such that every output symbol has at least this probability. */
    double minPrediction() const { return minPrediction_; }

    /**
     * Set the minimum probability for any observation symbol to the given `new_value`. The `prediction()` is normalized at each time-step such that every output symbol has at least this probability. */
    void minPrediction(double new_value) throw(std::invalid_argument) {
        if (new_value < 0) throw std::invalid_argument("minPrediction must be >= 0");
        else if (new_value > 1.0 / nO_) throw std::invalid_argument("minPrediction must be <= 1 / nOutputSymbols");
        minPrediction_ = new_value;
    }

    /**
     * Return the maximum allowed tolerance between the unnormalized and normalized prediction vector before invoking a `setBack()`. The tolerance is computed as 1.5 * `nOutputSymbols()` * squared norm of the difference. See also `normalizePrediction()`. */
    double normalizationTolerance() const { return normalizationTolerance_; }

    /**
     * Set the maximum allowed tolerance between the unnormalized and normalized prediction vector before invoking a `setBack()` to the given `new_value`. The tolerance is computed as 1.5 * `nOutputSymbols()` * squared norm of the difference. See also `normalizePrediction()`. */
    void normalizationTolerance(double new_value) throw(std::invalid_argument) {
        if (new_value < 0) throw std::invalid_argument("normalizationTolerance must be >= 0");
        normalizationTolerance_ = new_value;
    }

    /**
     * Return the maximum number of steps to "replay" during a `setBack()` operation. */
    int maxSetback() const { return maxSetback_; }

    /**
     * Set the maximum number of steps to "replay" during a `setBack()` operation to `new_value`. */
    void maxSetback(int new_value) throw(std::invalid_argument) {
        if (new_value < 0) throw std::invalid_argument("maxSetback must be >= 0");
        historyLength_ = std::min((long)(maxSetback_), historyLength_);
        maxSetback_ = (unsigned int)new_value;
        out_buf_.resize(maxSetback_, 0);
        in_buf_.resize(maxSetback_, 0);
    }

    /**
     * Return the smallest probability value considered as non-zero. If an observation symbol is encountered that this model predicts to have a probability lower than `impossibilityThreshold()`, this is considered as an impossible event and a `setBack()` must be performed, as the state `wt()` can no longer be normalized. */
    double impossibilityThreshold() const { return impossibilityThreshold_; }

    /**
     * Set the smallest probability value considered as non-zero to the given `new_value`. If an observation symbol is encountered that this model predicts to have a probability lower than `impossibilityThreshold()`, this is considered as an impossible event and a `setBack()` must be performed, as the state `wt()` can no longer be normalized. */
    void impossibilityThreshold(double new_value) throw(std::invalid_argument) {
        if (new_value < 0) throw std::invalid_argument("impossibilityThreshold must be >= 0");
        impossibilityThreshold_ = new_value;
    }
    //@} </editor-fold>

    /* <editor-fold desc="Basic OOM functionality"> */ /** @name Basic OOM functionality */ //@{
    /**
     * Return the (normalized) current state vector \f$\omega_t\f$ */
    const VectorXd &wt() const { return wt_; }

    /**
     * Return the most recent input-output history that is relevant for stabilization purposes. */
    Sequence history() const {
        int available_history_length = (int) (historyLength_ < maxSetback_ ? historyLength_ : maxSetback_);
        Sequence h(available_history_length, nOutputSymbols(), nInputSymbols());
        for (int i = 0; i < available_history_length; ++i) { h.o(available_history_length -1 -i, out_buf_[i]); }
        if (nU_ != 0) {
            for (int i = 0; i < available_history_length; ++i) { h.u(available_history_length -1 -i, in_buf_[i]); }
        }
        return h;
    }

    /**
     * Set the current state to the (normalized) given vector \c new_value, and (optionally) specify a given input-output `history` (relevant for stabilization). In the case of an output-only `Oom` this automatically calls `condition()`. */
    void wt(const VectorXd &new_value, Sequence history = Sequence()) {
        historyLength_ = history.length() < maxSetback_ ? history.length() : maxSetback_;
        for (int i = 0; i < historyLength_; ++i) {
            out_buf_[i] = history.o(-(i+1)); // uses -1 based negative indexing
            if (nU_ != 0) { in_buf_[i] = history.u(-(i+1)); }
        }
        didSetback_ = false;
        wt_ = new_value;
        double sig_wt = sig_ * wt_;
        if (sig_wt < impossibilityThreshold_) setBack();
        else wt_ /= sig_wt;
        if (nU_ == 0) condition();
    }

    /**
     * Reset the `Oom` to its initial state and `resetStabilizationStatistics()`. */
    void reset() {
        historyLength_ = 0;
        didSetback_ = false;
        wt_ = w0_;
        if (nU_ == 0) condition();
        nSetback_ = 0;
        nFixPrediction_ = 0;
        nImpossible_ = 0;
    }

    /**
     * Update the `Oom` state according to the input-output pair (\c u,\c o), normalize, and in the case of an output-only `Oom`, additionally call `condition()`.
     *
     * That is, first set `wt` to `tau(o,u)` * `wt()`, and then attempt to normalize to `wt()` / ( `sig()` * `wt()` ) if `sig()` * `wt()` is greater than the `impossibilityThreshold()`, else perform `setback()` operations. */
    void update(Symbol o, Symbol u = 0) {
        // update state
        temp_dim_.noalias() = tau_(o, u) * wt_;
        wt_ = temp_dim_;
        out_buf_.push_front(o);
        out_buf_.pop_back();
        ++historyLength_;
        if (nU_ != 0) {
            in_buf_.push_front(u);
            in_buf_.pop_back();
        }
        didSetback_ = false;
        // normalize state
        double sig_wt = sig_ * wt_;
        if (sig_wt < impossibilityThreshold_) setBack();
        else wt_ /= sig_wt;
        // if not IO, condition, i.e., get and fix probability vector:
        if (nU_ == 0) condition();
    }

    /**
     * Compute and normalize the prediction vector of the next output symbol probabilities according to the current state `wt()` and input `u`, i.e., compute P( · | `u`, `wt()` ).
     *
     * Note that this function calls `setback()` if the `normalizationTolerance()` is exceeded. Therefore, calling e.g., `condition(0)` after `condition(1)` may not give the same result as calling just `condition(0)`. */
    void condition(Symbol u = 0) {
        P_.noalias() = sig_tau_(u) * wt_;
        double fP = 0;
        while (((fP = normalizePrediction()) > normalizationTolerance_) and setBack())
            P_.noalias() = sig_tau_(u) * wt_;
        if (fP > fixPredictionMargin_) ++nFixPrediction_;
    }

    /**
     * Return the current prediction vector of the next output symbol probabilities. In the case of an input-output `Oom` these probabilities depend on the current input symbol u_t, so `condition(u_t)` *must* have been called first. */
    const VectorXd &prediction() const { return P_; }

    SWIGCODE(%feature ("kwargs") sample;)
    /**
     * Sample, using the given `randomSource`, a sequence of given `length` from the `Oom`, in the case of inputs together with the given input `policy` (by default iud inputs), starting from the **current** state `wt()`.
     *
     * For each time-step an observation is sampled from the `prediction()` vector raised element-wise to the power `exponent` (default 1). Higher `exponent` values add a bias towards the most likely sequences, while lower values bias towards uniformly distributed sequences. */
    Sequence sample(unsigned long length, const Random &randomSource = Random(), const Policy &policy = Policy(), double exponent = 1) {
        if (policy.nU_ > nU_) { throw std::invalid_argument("incompatible policy"); }
        Sequence seq(length, nO_, nU_);
        Symbol u = 0, o = 0;
        for (unsigned long l = 0; l < length; l++) {
            LOOP_PROGRESS("Generating sequence", l, length)
            if (nU_ != 0) {
                u = policy.nU_ == 0 ? randomSource.integer(nInputSymbols()) : policy.u(wt(), randomSource);
                condition(u);
            }
            if (exponent != 1) { P_ = P_.array().pow(exponent); }
            o = randomSource.sample(P_);
            if (nU_ != 0) { seq.u(l, u); }
            seq.o(l, o);
            update(o, u);
        }
        LOOP_DONE("Generating sequence")
        return seq;
    }

    /**
     * If `reset` is `true` (default), perform a state `reset()` first\. Then return the value for the given output symbol `z` of the prediction function for the state `wt()`, i.e., the "probability" P( `z` | `wt()` ), and update the state. */
    double f(Symbol z, bool reset = true) throw(std::invalid_argument) {
        if (nU_ != 0) throw std::invalid_argument("input symbol required");
        return f(z, 0, reset);
    }

    /**
     * If `reset` is `true` (default), perform a state `reset()` first\. Then return the value for the given input-output pair (`u`,`o`) of the prediction function for the state `wt()`, i.e., the "probability" P( `o` | `u`, `wt()` ), and update the state. In the case of an output-only `Oom`, the input `u` is simply ignored. */
    double f(Symbol o, Symbol u, bool reset = true) {
        if (reset) { this->reset(); }
        if (nU_ != 0) condition(u);
        double val = prediction()(o);
        update(o, u);
        return val;
    }

    /**
     * If `reset` is `true` (default), perform a state `reset()` first\. Then return the value for the given `sequence` of the prediction function for the state `wt()`, i.e., the "probability" P( `sequence` | `wt()` ), and update the state. */
    double f(const Sequence &sequence, bool reset = true) {
        if (reset) { this->reset(); }
        double val = 1;
        for (unsigned long pos = 0; pos < sequence.length(); ++pos)
            val *= f(sequence.o(pos), ((nU_ != 0) ? sequence.u(pos) : 0), false);
        return val;
    }

    /**
     * If `reset` is `true` (default), perform a state `reset()` first\. Then return the matrix of prediction function values \f$[ f(x y) ]_{y \in Y, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences. */
    MatrixXd f(const Sequences &Y, const Sequences &X, bool reset = true) {
        return f(Y, Sequence(0, nO_, nU_), X, reset);
    }

    /**
     * If `reset` is `true` (default), perform a state `reset()` first\. Then return the matrix of prediction function values \f$[ f(x z y) ]_{y \in Y, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences for a given output symbol `z`. */
    MatrixXd f(const Sequences &Y, Symbol z, const Sequences &X, bool reset = true) {
        if (!isIO()) return f(Y, Sequence(std::vector<Symbol>{z}, nO_, nU_), X, reset);
        else         return MatrixXd();
    }

    /**
     * If `reset` is `true` (default), perform a state `reset()` first\. Then return the matrix of prediction function values \f$[ f(x z y) ]_{y \in Y, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences for a given input-output symbol pair z = (`u`, `o`). In the case of an output-only `Oom`, the input `u` is simply ignored. */
    MatrixXd f(const Sequences &Y, Symbol o, Symbol u, const Sequences &X, bool reset = true) {
        if (isIO()) return f(Y, Sequence(std::vector<Symbol>{u, o}, nO_, nU_), X, reset);
        else        return f(Y, Sequence(std::vector<Symbol>{o}, nO_, nU_), X, reset);
    }

    /**
     * If `reset` is `true` (default), perform a state `reset()` first\. Then return the matrix of prediction function values \f$[ f(x s y) ]_{y \in Y, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences for a given `Sequence` `s`. */
    MatrixXd f(const Sequences &Y, Sequence s, const Sequences &X, bool reset = true) {
        if (reset) { this->reset(); }
        double f_temp;
        VectorXd wt_temp;
        Sequence history_temp;
        VectorXd wt_old = wt_;
        Sequence history_old = history();
        unsigned long Frows = Y.size(), Fcols = X.size();
        MatrixXd F = MatrixXd::Zero(Frows, Fcols);
        for (unsigned int j = 0; j < Fcols; ++j) {
            wt(wt_old, history_old);
            f_temp = f(X[j], false);
            wt_temp = wt_;
            history_temp = history();
            for (unsigned int i = 0; i < Frows; ++i) {
                wt(wt_temp, history_temp);
                f_temp *= f(s, false);
                F.coeffRef(i, j) = f_temp * f(Y[i], false);
            }
        }
        return F;
    }


    SWIGCODE(%feature ("kwargs") log2_f;)
    /**
     * If `reset` is `true` (default), perform a state `reset()` first\. Then return the log_2 of the prediction function for the given `sequence` given the state `wt()`, i.e., log_2 f( `sequence` | `wt()` ), and update the state.
     *
     * To deal gracefully with observations that have a prediction value below `impossibilityThreshold()` at some time step but occur nevertheless in the `sequence`, the prediction for this occurrence is treated as having a probability of `impossibilityThreshold()`. Every time this happens, the counter  `nImpossible_` is incremented. Note that this problem can be avoided by increasing `minPrediction()` above `impossibilityThreshold()`. */
    double log2_f(const Sequence &sequence, bool reset = true) {
        if (reset) { this->reset(); }
        double LL = 0;
        double val;
        for (unsigned long pos = 0; pos < sequence.length(); ++pos) {
            val = f(sequence.o(pos), sequence.u(pos), false);
            if (val < impossibilityThreshold_) {
                LL += log2(impossibilityThreshold_);
                nImpossible_++;
            } else {
                LL += log2(val);
            }
        }
        return LL;
    }

    /**
     * Return the log2-likelihood of the `Oom` for the given `sequence`. That is, return -`log2_f(sequence, true)` / `sequence.length()`. */
    double l2l(const Sequence &sequence) { return -log2_f(sequence, true) / sequence.length(); }

    /**
     * Return entropy of this `Oom` estimated on a sample sequence of the given `sample_length`.
    */
    double entropy(unsigned long sample_length, const Random &randomSource = Random(), const Policy &policy = Policy()) {
        if (policy.nU_ > nU_) { throw std::invalid_argument("incompatible policy"); }
        Symbol u = 0, o = 0;
        double entropy_estimate = 0;
        for (unsigned long l = 0; l < sample_length; l++) {
            LOOP_PROGRESS("Estimating entropy", l, sample_length)
            if (nU_ != 0) {
                u = policy.nU_ == 0 ? randomSource.integer(nInputSymbols()) : policy.u(wt(), randomSource);
                condition(u);
            }

            entropy_estimate -= (Eigen::log((P_.array() > 0).select(P_.array(), 1.0)) * P_.array()).sum();
            o = randomSource.sample(P_);
            update(o, u);
        }
        LOOP_DONE("Estimating entropy")
        return entropy_estimate / log(2.0) / sample_length;
    }

    /**
     * Return the cross-entropy of the best `k`-order Markov model approximation estimated on the given sample `sequence`.
    */
    double crossEntropyOfKOrderMarkovApproximation(int k, const Sequence& sequence);


    SWIGCODE(%feature ("kwargs") averageOneStepPredictionError;)
    /**
     * Return the average one-step squared prediction error computed along the given sample `sequence` according to a correct `Oom` `trueModel`, after first performing a state `reset()` operation on both this `Oom` and the `trueModel`. Their states are updated accordingly. */
    double averageOneStepPredictionError(const Sequence &sequence, Oom &trueModel) {
        reset();
        trueModel.reset();
        double sqr_err = 0;
        for (unsigned long pos = 0; pos < sequence.length(); pos++) {
            Symbol u = sequence.u(pos);
            Symbol o = sequence.o(pos);
            condition(u);
            trueModel.condition(u);
            sqr_err += (prediction() - trueModel.prediction()).squaredNorm() / nO_;
            update(o, u);
            trueModel.update(o, u);
        }
        return sqr_err / sequence.length();
    }

    /**
     * Return the `float` matrix of the `sequence.length()` states (in its columns) occurring during the computation of `f(sequence, reset)`. */
    MatrixXf harvestStates(const Sequence &sequence, bool reset = true) {
        if (reset) { this->reset(); }
        long sequence_length = sequence.length();
        Symbol u = 0;
        MatrixXf W(dim_, sequence_length + 1);
        W.col(0) = wt().cast<float>();
        for (long i = 0; i < sequence_length; ++i) {
            if (nU_ != 0) {
                u = sequence.u(i);
                condition(u);
            }
            update(sequence.o(i), u);
            W.col(i + 1) = wt().cast<float>();
        }
        return W;
    }
    //@} </editor-fold>

    /* <editor-fold desc="Transformation functions"> */ /** @name Transformation functions */ //@{

    // TODO: review `reverse()`, `transform()` and `conjugate()`

    SWIGCODE(%feature ("kwargs") stationaryState;)
    /**
     * Return the stationary state (in the case of an input-output `Oom` according to the given `policy`), computed by the power method with at most `maxIterations` number of iterations. */
    Eigen::VectorXd stationaryState(const Policy& policy = Policy(), int maxIterations = 10000) const {
        if (policy.nU_ > nU_) { throw std::invalid_argument("incompatible policy"); }
        Array<MatrixXd, Dynamic, 1> t(nU_ == 0 ? 1 : nU_);
        for (Symbol u = 0; u < (nU_ == 0 ? 1 : nU_); ++u) {
            t(u) = MatrixXd::Zero(dim_, dim_);
            for (Symbol o = 0; o < nO_; ++o) {
                t(u) += tau(o, u);
            }
        }
        MatrixXd tt;
        VectorXd w = w0();
        VectorXd w_temp = w;
        VectorXd inputProbs = VectorXd::Ones(nU_ == 0 ? 1 : nU_) / (nU_ == 0 ? 1 : nU_);
        double error = 1;
        double lambda = 1;
        int iteration = 0;
        while (error > fabs(lambda) * std::numeric_limits<double>::epsilon()
               and iteration < maxIterations) {
            iteration++;
            w = w_temp;
            w.normalize();
            if (policy.nU_ != 0) { inputProbs = policy.p(w); }
            tt.setZero(dim_, dim_);
            for (Symbol u = 0; u < (nU_ == 0 ? 1 : nU_); ++u) {
                tt += inputProbs(u) * t(u);
            }
            w_temp.noalias() = tt * w;
            lambda = (double)(w.transpose() * w_temp);
            error = (w_temp - lambda * w).norm();
        }
        w /= (double)(sig_ * w);
        return w;
    }

    /**
     * Return the "reverse" of this `Oom`. */
    std::shared_ptr<Oom> reverse(bool normalize = true) const {
        auto room = std::make_shared<Oom>(*this); // copy of oom;
        MatrixXd rho = MatrixXd::Identity(dim_, dim_);
        MatrixXd rhoInv = MatrixXd::Identity(dim_, dim_);
        if (normalize) {
            rho = w0().asDiagonal().inverse(); // rho * w0 = (1,...,1)
            rho /= dimension(); // it should also work without this
            rhoInv = dimension() * w0().asDiagonal();
        }
        for (Symbol a = 0; a < (nU_ == 0 ? 1 : nU_); ++a)
            for (Symbol o = 0; o < nOutputSymbols(); ++o)
                room->tau(o, a, (rho * tau(o, a) * rhoInv).transpose());
        room->w0((sig() * rhoInv).transpose());
        room->sig((rho * w0()).transpose());
        room->initialize();
        return room;
    }

    /**
     * Transform this `Oom` to an equivalent `Oom` that has given `sig` and `w0` as parameters for `sig()` and `w0()`\. This will only yield an (equivalent) `Oom` if `sig` * `w0` = 1. */
    void transform(const RowVectorXd &sig, const VectorXd &w0 = VectorXd::Zero(0)) {
        VectorXd e1 = VectorXd::Zero(dimension());
        e1(0) = 1;
        VectorXd v = sig_.transpose() / sig_.norm() - e1;
        MatrixXd H = MatrixXd::Identity(dimension(), dimension());
        if (not v.isZero()) {
            v /= v.norm();
            H -= 2 * v * v.transpose();
        }
        H.col(0) = w0_;
        v = sig.transpose() / sig.norm() - e1;
        MatrixXd H2 = MatrixXd::Identity(dimension(), dimension());
        if (not v.isZero()) {
            v /= v.norm();
            H2 -= 2 * v * v.transpose();
        }
        H2.col(0) = w0;
        conjugate(H2 * H.inverse(), H * H2.inverse());
    }

    /**
     * Conjugate this `Oom` by the given matrices `rho` and `rhoInv`.
     *
     * That is, set:
     * - `w0()` = `rho` * `w0()`
     * - `tau(o,u)` = `rho` * `tau(o,u)` * `rhoInv`
     * - `sig()` = `sig()` * `rhoInv`.
     */
    void conjugate(const MatrixXd &rho, const MatrixXd &rhoInv) {
        w0_ = rho * w0_;
        for (Symbol u = 0; u < (nU_ == 0 ? 1 : nU_); u++)
            for (Symbol o = 0; o < nO_; o++)
                tau_(o, u) = rho * tau_(o, u) * rhoInv;
        sig_ = sig_ * rhoInv;
        dim_ = (int) rho.rows();
        initialize();
    }
    //@} </editor-fold>

    bool operator==(const Oom &other) const {
        if (dimension() != other.dimension() or nOutputSymbols() != other.nOutputSymbols() or nInputSymbols() != other.nInputSymbols()) return false;
        if (sig() != other.sig() or w0() != other.w0()) return false;
        for (Symbol o = 0; o < nOutputSymbols(); ++o)
            for (Symbol u = 0; u < ((nInputSymbols() == 0) ? 1 : nInputSymbols()); ++u)
                if (tau(o, u) != other.tau(o, u)) return false;
        if (maxSetback() != other.maxSetback()
            or minPrediction_ != other.minPrediction_
            or normalizationTolerance_ != other.normalizationTolerance_
            or impossibilityThreshold_ != other.impossibilityThreshold_) return false;
        return true;
    }

    bool operator!=(const Oom &other) const { return not operator==(other); }

    /* <editor-fold desc="IO-functions"> */ /** @name IO-functions */ //@{
    INSERT_JSON_IO_FUNCTIONS()

    /** return a representation to display in interactive python. */
    std::string repr() const {
        std::stringstream os;
        os << "OOM" << " nU: " << nU_ << " nO: " << nO_ << " dim: " << dim_;
        return os.str();
    }
    //@} </editor-fold>

    /* <editor-fold desc="Internals"> */ /** @name Internalals\. Use only if you know what you are doing! */ //@{
    /**
     * Attempt to fix the prediction vector of the next output symbol probabilities \f$P(\cdot|u_t, \omega_t)\f$ such that all probabilities are at least \c minPrediction_ and the probabilities sum to one.\ Return a measure of the required change to the prediction vector: 1.5 * nO() * squared norm of the difference. */
    double normalizePrediction() {
        temp_nO_ = P_;
        P_.array() -= minPrediction_;
        P_ = P_.cwiseMax(0);
        double P_sum = P_.sum();
        if (P_sum < impossibilityThreshold_) P_.setConstant(1. / nO_);
        else P_ = P_.array() * (1 - nO_ * minPrediction_) / P_sum + minPrediction_;
        return sqrt((P_ - temp_nO_).squaredNorm() / nO_);
        // This used to be (P_ - oldP_).squaredNorm() * 1.5 * nO_;
    }

    /**
     * Attempt to perform a state setback operation for at most \c maxSetback_ time-steps.\ Note that calling this method repeatedly will attempt a setback for a shorter history each time.\ Return \c true if a setback could be performed. */
    bool setBack() {
        if (historyLength_ == 0) {
            wt_ = w0_;
            return false;
        } // nothing to setback!
        if (!didSetback_) {
            nSetback_++;
            didSetback_ = true;
        }
        double sig_wt = 0;
        historyLength_ = std::min((long) (maxSetback_ + 1), historyLength_); // --historylength next...
        while ((sig_wt < impossibilityThreshold_) and (--historyLength_ > 0)) {
            // we decrease history by one, otherwise we just get the same state back
            wt_ = w0_;
            for (int i = 0; i < historyLength_; ++i) {
                // replay the history
                temp_dim_ = tau_(out_buf_[historyLength_ - 1 - i], in_buf_[historyLength_ - 1 - i]) * wt_;
                wt_ = temp_dim_;
                sig_wt = sig_ * wt_;
                if (sig_wt < impossibilityThreshold_) break; // still bad states: try with shorter history...
                wt_ /= sig_wt;
            }
        }
        if (historyLength_ == 0) wt_ = w0_;
        return true;
    }
    //@} </editor-fold>

    double fixPredictionMargin_ = 0.01; /**< the largest tolerated error of the prediction vector before the normalization is considered a fixing event */
    int nSetback_ = 0; /**< the number of times that the OOM state needed to be fixed by a setback operation (introduces a considerable error) */
    int nFixPrediction_ = 0; /**< the number of times that the predicted probabilities were adjusted by more than \c fixPredictionMargin_. */
    int nImpossible_ = 0; /**< the number of times that a symbol was encountered that has probability smaller than \c impossibleProbMargin_ according to the OOM */

private:
    int dim_; /**< the model dimension */
    int nO_; /**< the size of the output alphabet */
    int nU_; /**< the size of the input alphabet */
    RowVectorXd sig_; /**< the probability functional vector \f$\sigma\f$ */
    Array<MatrixXd, Dynamic, Dynamic> tau_; /**< the array of observable operators (of size \c nO x \c nU). The operator corresponding to the input-observation pair (u,o) is given by tau(o,u). (Read: "the observable operator for the observation o given the input u"). This is convenient when dealing with standard OOMs (\c nU == 0), since one can address the operators as tau(o) in this case, as the action u defaults to 0). */
    VectorXd w0_; /**< the initial state */
    VectorXd wt_; /**< the current state (at time t) */
    double minPrediction_ = 0 /*1e-5*/; /**< the smallest probability to assign to any next-symbol prediction when normalizing predictions */
    double normalizationTolerance_ = 0.3; /**< the largest tolerated error of the prediction vector before a state setback is performed */
    double impossibilityThreshold_ = 1e-12; /**< the smallest number reasonably considered non-zero (e.g., for division purposes) */
    unsigned int maxSetback_ = 0; /**< the maximum number of steps to "replay" during a setBack operation */

    bool didSetback_ = false; /**< \c true if a setback operation has been performed in the current time step */
    Array<MatrixXd, 1, Dynamic> sig_tau_;
    std::deque<Symbol> in_buf_, out_buf_;
    VectorXd P_; /**< the conditional next-symbol probabilities given the current input: \f$P(o)=P(o|u_t, w_t)=\sigma*tau(o,u)*w_t\f$. */
    long historyLength_ = 0;
    VectorXd temp_dim_, temp_nO_; /**< preallocated vectors of the named size to use as temporaries */

    #ifndef SWIG

    template<class Archive>
    void save(Archive &ar) const {
        const std::string type = "OOM";
        CEREALIZE(ar, type, Type);
        CEREALIZE(ar, nU_, nU);
        CEREALIZE(ar, nO_, nO);
        CEREALIZE(ar, dim_, dim);
        CEREALIZE(ar, sig_, sig);
        CEREALIZE(ar, tau_, tau);
        CEREALIZE(ar, w0_, w0);
        CEREALIZE(ar, minPrediction_, minPrediction);
        CEREALIZE(ar, normalizationTolerance_, normalizationTolerance);
        CEREALIZE(ar, maxSetback_, maxSetback);
        CEREALIZE(ar, impossibilityThreshold_, impossibilityThreshold);
    }

    template<class Archive>
    void load(Archive &ar) {
        std::string type;
        CEREALIZE(ar, type, Type);
        CEREALIZE(ar, nU_, nU);
        CEREALIZE(ar, nO_, nO);
        CEREALIZE(ar, dim_, dim);
        CEREALIZE(ar, sig_, sig);
        CEREALIZE(ar, tau_, tau);
        CEREALIZE(ar, w0_, w0);
        stabilization(-1, -1, -1, -1, "none");
        CEREALIZE_OPTIONAL(ar, minPrediction_, minPrediction);
        /* For backwards compatibility: */
        try { CEREALIZE(ar, normalizationTolerance_, normalizationTolerance); }
        catch (...) { CEREALIZE_OPTIONAL(ar, normalizationTolerance_, maxPredictionError); }
        CEREALIZE_OPTIONAL(ar, maxSetback_, maxSetback);
        CEREALIZE_OPTIONAL(ar, impossibilityThreshold_, impossibilityThreshold);
        initialize();
    }

    #endif // SWIG
};

} // namespace tom
