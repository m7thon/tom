namespace tom {

SWIGCODE(%feature("python:slot", "tp_repr", functype="reprfunc") Estimator::repr;)

/**
 * This class computes estimates for \f$f( x )\f$ and corresponding variance estimates for sequences \f$x\f$ based on a suffix tree representation of a sample sequence.
 */
class Estimator {
public:
    /**
     * Create an `Estimator` for a sample sequence data given by a suffix tree representation `stree`. */
    Estimator(const std::shared_ptr<stree::STree> &stree)
            : state_(stree), stree_(stree) {
        nO_ = stree_->sequence().nOutputSymbols();
        nU_ = stree_->sequence().nInputSymbols();
        len_ = (stree_->sequence().length());
        uProbs_ = VectorXd::Ones(std::max(1, nU_));
        for (Symbol u = 0; u < nU_; ++u) {
            state_.pos_.setRoot();
            state_.pos_.toSymbol(u);
            uProbs_(u) = (double) (state_.pos_.count()) / len_;
        }
        state_.pos_.setRoot();
        regularization(-1, -1, -1, -1, "default");
    }

    /** Return the data sequence. */
    Sequence sequence() const {
        return stree_->sequence();
    }

    //<editor-fold desc="Estimates by f, v and fv">
    /**
     * Return an estimate of f( `z` ) for the given output symbol `z`. */
    double f(Symbol z) throw(std::invalid_argument) {
        if (nU_ != 0) throw std::invalid_argument("input symbol required");
        return f(Sequence(std::vector<Symbol>{z}, nO_, nU_));
    }

    /**
     * Return an estimate of f( z ) for the given input-output symbol pair z = (`u`, `o`). In the case of an output-only system, the input `u` is simply ignored. */
    double f(Symbol o, Symbol u) {
        if (nU_ == 0) return f(Sequence(std::vector<Symbol>{o}, nO_, nU_));
        else          return f(Sequence(std::vector<Symbol>{u, o}, nO_, nU_));
    }

    /**
     * Return an estimate of f( `sequence` ). */
    double f(const Sequence &sequence) {
        estimateVariance_ = false;
        state_.reset();
        extendBy(sequence);
        eval();
        return state_.f_;
    }

    /**
     * Return the matrix of estimates for \f$[ f( x y ) ]_{y \in X, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences. */
    MatrixXd f(const Sequences &Y, const Sequences &X) { return f(Y, Sequence(0, nO_, nU_), X); }

    /**
     * Return the matrix of estimates for \f$[ f( x z y ) ]_{y \in X, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences for a given output symbol `z`. */
    MatrixXd f(const Sequences &Y, Symbol z, const Sequences &X) {
        if (nU_ == 0) return f(Y, Sequence(std::vector<Symbol>{z}, nO_, nU_), X);
        return MatrixXd();
    }

    /**
     * Return the matrix of estimates for \f$[ f( x z y ) ]_{y \in X, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences for a given input-output symbol pair z = (`u`, `o`). In the case of an output-only system, the input `u` is simply ignored. */
    MatrixXd f(const Sequences &Y, Symbol o, Symbol u, const Sequences &X) {
        if (nU_ != 0) return f(Y, Sequence(std::vector<Symbol>{u, o}, nO_, nU_), X);
        else          return f(Y, Sequence(std::vector<Symbol>{o}, nO_, nU_), X);
    }

    /**
     * Return the matrix of estimates for \f$[ f( x s y ) ]_{y \in X, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences for a given `Sequence` `s`. */
    MatrixXd f(const Sequences &Y, const Sequence s, const Sequences &X) {
        unsigned long rows = Y.size(), cols = X.size();
        MatrixXd F(rows, cols);
        State state(stree_);
        estimateVariance_ = false;
        for (unsigned long j = 0; j < cols; ++j) {
            state_.reset();
            extendBy(X[j]);
            state = state_;
            for (unsigned long i = 0; i < rows; ++i) {
                state_ = state;
                extendBy(s);
                extendBy(Y[i]);
                eval();
                F.coeffRef(i, j) = state_.f_;
            }
        }
        return F;
    }

    /**
     * Return a variance estimate for the estimate of f( `z` ) for the given output symbol `z`. */
    double v(Symbol z) throw(std::invalid_argument) {
        if (nU_ != 0) throw std::invalid_argument("input symbol required");
        return v(Sequence(std::vector<Symbol>{z}, nO_, nU_));
    }

    /**
     * Return a variance estimate for the estimate of f( z ) for the given input-output symbol pair z = (`u`, `o`). In the case of an output-only system, the input `u` is simply ignored. */
    double v(Symbol o, Symbol u) {
        if (nU_ != 0) return v(Sequence(std::vector<Symbol>{u, o}, nO_, nU_));
        else          return v(Sequence(std::vector<Symbol>{o}, nO_, nU_));
    }

    /**
     * Return a variance estimate for the estimate of f( `sequence` ). */
    double v(const Sequence &sequence) {
        estimateVariance_ = true;
        state_.reset();
        extendBy(sequence);
        eval();
        return state_.v_;
    }

    /**
     * Return the matrix of element-wise variance estimates corresponding to the estimates for \f$[ f( x y ) ]_{y \in X, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences. */
    MatrixXd v(const Sequences &Y, const Sequences &X) { return v(Y, Sequence(0, nO_, nU_), X); }

    /**
     * Return the matrix of element-wise variance estimates corresponding to the estimates for \f$[ f( x z y ) ]_{y \in X, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences for a given output symbol `z`. */
    MatrixXd v(const Sequences &Y, Symbol z, const Sequences &X) {
        if (nU_ == 0) return v(Y, Sequence(std::vector<Symbol>{z}, nO_, nU_), X);
        return MatrixXd();
    }

    /**
     * Return the matrix of element-wise variance estimates corresponding to the estimates for \f$[ f( x z y ) ]_{y \in X, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences for a given input-output symbol pair z = (`u`, `o`). In the case of an output-only system, the input `u` is simply ignored. */
    MatrixXd v(const Sequences &Y, Symbol o, Symbol u, const Sequences &X) {
        if (nU_ != 0) return v(Y, Sequence(std::vector<Symbol>{u, o}, nO_, nU_), X);
        else          return v(Y, Sequence(std::vector<Symbol>{o}, nO_, nU_), X);
    }

    /**
     * Return the matrix of element-wise variance estimates corresponding to the estimates for \f$[ f( x s y ) ]_{y \in X, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences for a given `Sequence` `s`. */
    MatrixXd v(const Sequences &Y, const Sequence s, const Sequences &X) {
        MatrixXd F, V;
        fv(F, V, Y, s, X);
        return V;
    }

    SWIGCODE(%apply MatrixXd& OUTPUT { MatrixXd& F, MatrixXd& V };)
    SWIGCODE(%apply double& OUTPUT { double& f, double& v };)

    /**
     * Return
     * \if PY
     * in a tuple (`f`, `v`)
     * \else
     * (in the output arguments `f` and `v`)
     * \endif
     * an estimate `f` of f( `z` ) for the given output symbol `z` together with the corresponding variance estimate `v`. */
    C1(void) PY2(tuple<double, double>)
    fv(C3(double &f, double &v,) Symbol z) throw(std::invalid_argument) {
        if (nU_ != 0) throw std::invalid_argument("input symbol required");
        fv(f, v, Sequence(std::vector<Symbol>{z}, nO_, nU_));
    }

    /**
     * Return
     * \if PY
     * in a tuple (`f`, `v`)
     * \else
     * (in the output arguments `f` and `v`)
     * \endif
     * an estimate `f` of f( z ) for the given input-output symbol pair z = (`u`, `o`) together with the corresponding variance estimate `v`. In the case of an output-only system, the input `u` is simply ignored. */
    C1(void) PY2(tuple<double, double>)
    fv(C3(double &f, double &v,) Symbol o, Symbol u) {
        if (nU_ == 0) fv(f, v, Sequence(std::vector<Symbol>{o}, nO_, nU_));
        else          fv(f, v, Sequence(std::vector<Symbol>{u, o}, nO_, nU_));
    }

    /**
     * Return
     * \if PY
     * in a tuple (`f`, `v`)
     * \else
     * (in the output arguments `f` and `v`)
     * \endif
     * an estimate of f( `sequence` ) together with the corresponding variance estimate `v`. */
    C1(void) PY2(tuple<double, double>)
    fv(C3(double &f, double &v,) const Sequence &sequence) {
        estimateVariance_ = true;
        state_.reset();
        extendBy(sequence);
        eval();
        f = state_.f_;
        v = state_.v_;
    }

    /**
     * Return
     * \if PY
     * in a tuple (`F`, `V`)
     * \else
     * (in the output arguments `F` and `V`)
     * \endif
     * the matrix `F` of estimates for \f$[ f( x y ) ]_{y \in X, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences, together with the corresponding matrix `V` of element-wise variance estimates for the estimates returned in `F`. */
    C1(void) PY2(tuple<MatrixXd, MatrixXd>)
    fv(C3(MatrixXd &F, MatrixXd &V,) const Sequences &Y, const Sequences &X) {
        fv(F, V, Y, Sequence(0, nO_, nU_), X);
    }

    /**
     * Return
     * \if PY
     * in a tuple (`F`, `V`)
     * \else
     * (in the output arguments `F` and `V`)
     * \endif
     * the matrix `F` of estimates for \f$[ f( x z y ) ]_{y \in X, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences for a given output symbol `z`, together with the corresponding matrix `V` of element-wise variance estimates for the estimates returned in `F`. */
    C1(void) PY2(tuple<MatrixXd, MatrixXd>)
    fv(C3(MatrixXd &F, MatrixXd &V,) const Sequences &Y, Symbol z, const Sequences &X) {
        if (nU_ == 0) fv(F, V, Y, Sequence(std::vector<Symbol>{z}, nO_, nU_), X);
    }

    /**
     * Return
     * \if PY
     * in a tuple (`F`, `V`)
     * \else
     * (in the output arguments `F` and `V`)
     * \endif
     * the matrix `F` of estimates for \f$[ f( x z y ) ]_{y \in X, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences for a given input-output symbol pair z = (`u`, `o`), together with the corresponding matrix `V` of element-wise variance estimates for the estimates returned in `F`. In the case of an output-only system, the input `u` is simply ignored. */
    C1(void) PY2(tuple<MatrixXd, MatrixXd>)
    fv(C3(MatrixXd &F, MatrixXd &V,) const Sequences &Y, Symbol o, Symbol u, const Sequences &X) {
        if (nU_ == 0) fv(F, V, Y, Sequence(std::vector<Symbol>{o}, nO_, nU_), X);
        else          fv(F, V, Y, Sequence(std::vector<Symbol>{u, o}, nO_, nU_), X);
    }

    /**
     * Return
     * \if PY
     * in a tuple (`F`, `V`)
     * \else
     * (in the output arguments `F` and `V`)
     * \endif
     * the matrix `F` of estimates for \f$[ f( x s y ) ]_{y \in X, x \in X}\f$ with rows indexed by the given set `Y` of characteristic sequences and columns indexed by the given set `X` of indicative sequences for a given `Sequence` `s`, together with the corresponding matrix `V` of element-wise variance estimates for the estimates returned in `F`. */
    C1(void) PY2(tuple<MatrixXd, MatrixXd>)
    fv(C3(MatrixXd &F, MatrixXd &V,) const Sequences &Y, const Sequence s, const Sequences &X) {
        unsigned long rows = Y.size(), cols = X.size();
        F.resize(rows, cols);
        V.resize(rows, cols);
        State state(stree_);
        estimateVariance_ = true;
        for (unsigned long j = 0; j < cols; ++j) {
            state_.reset();
            extendBy(X[j]);
            state = state_;
            for (unsigned long i = 0; i < rows; ++i) {
                state_ = state;
                extendBy(s);
                extendBy(Y[i]);
                eval();
                F.coeffRef(i, j) = state_.f_;
                V.coeffRef(i, j) = state_.v_;
            }
        }
    }

    SWIGCODE(%clear MatrixXd& F;)
    SWIGCODE(%clear MatrixXd& V;)
    SWIGCODE(%clear double& f;)
    SWIGCODE(%clear double& v;)
    //</editor-fold>

    //<editor-fold desc="Regularization parameters">
    SWIGCODE(%feature ("kwargs") regularization;)
    SWIGCODE(%apply double *OUTPUT { double *nPseudoCountsOut, double *zConfidenceIntervalSizeOut, double *minimumVarianceOut, double *exponentOut };)
    /**
     Set (optional) and then return the regularization parameters for the `Estimator`. This is done as follows:
     - First, if a `preset` ("none" / "default") is specified, all regularization parameters are set accordingly.
     - Next, any non-default argument causes the corresponding regularization parameter to be set to the given value, while any argument left at its default value (-1) has no effect.
     - Finally, the current regularization parameters are returned
     \if PY
     in a tuple in the same order as they appear as function arguments. This allows writing python code such as:
     \verbatim
     old_params = estimator.regularization()
     ...
     estimator.regularization(*old_params)
     \endverbatim
     \else
     in the corresponding output arguments `...Out` (if not NULL).
     \endif
     */
    C1(void) PY1(tuple)
    regularization(C5(double *nPseudoCountsOut, double *zConfidenceIntervalSizeOut, double *minimumVarianceOut, double *exponentOut,)
                   double nPseudoCounts = -1,
                   double zConfidenceIntervalSize = -1,
                   double minimumVariance = -1,
                   double exponent = -1,
                   std::string preset = "") throw(std::invalid_argument) {
        if (preset == "default") {
            nPseudoCounts_           = 1;
            zConfidenceIntervalSize_ = 1;
            minimumVariance_         = 1;
            exponent_                = 0.5;
        } else if (preset == "none") {
            nPseudoCounts_           = 0;
            zConfidenceIntervalSize_ = 0;
            minimumVariance_         = 0;
            exponent_                = 1;
        } else if (preset != "") {
            throw std::invalid_argument("unrecognized preset");
        }
        if (nPseudoCounts           != -1) this->nPseudoCounts_ = nPseudoCounts;
        if (zConfidenceIntervalSize != -1) this->zConfidenceIntervalSize_ = zConfidenceIntervalSize;
        if (minimumVariance         != -1) this->minimumVariance_ = minimumVariance;
        if (exponent                != -1) this->exponent_ = exponent;
        if (nPseudoCountsOut) *nPseudoCountsOut = this->nPseudoCounts_;
        if (zConfidenceIntervalSizeOut) *zConfidenceIntervalSizeOut = this->zConfidenceIntervalSize_;
        if (minimumVarianceOut) *minimumVarianceOut = this->minimumVariance_;
        if (exponentOut) *exponentOut = this->exponent_;
    }

    SWIGCODE(%clear double *nPseudoCountsOut, double *zConfidenceIntervalSizeOut, double *minimumVarianceOut, double *exponentOut;)

#ifndef SWIG
#ifndef PY
    /**
     * Set the regularization parameters for the `Estimator`. This is a convenience function that simply calls `regularization(NULL, NULL, NULL, NULL, nPseudoCounts, zConfidenceIntervalSize, minimumVariance, exponent, preset);` */
    void regularization(double nPseudoCounts = -1,
                        double zConfidenceIntervalSize = -1,
                        double minimumVariance = -1,
                        double exponent = -1,
                        std::string preset = "") throw(std::invalid_argument) {
        regularization(NULL, NULL, NULL, NULL, nPseudoCounts, zConfidenceIntervalSize, minimumVariance, exponent, preset);
    }
#endif
#endif

    double nPseudoCounts() const { return nPseudoCounts_; }

    double zConfidenceIntervalSize() const { return zConfidenceIntervalSize_; }

    double minimumVariance() const { return minimumVariance_; }

    double exponent() const { return exponent_; }
    //</editor-fold>

    VectorXd uProbs_; ///< the input Symbol probabilities for the case of an iid ("blind") input policy.\ These are estimated on constructing the \c Estimator and may be overwritten if they are known exactly.

private:
    int nO_;              ///< the size of the output alphabet
    int nU_;              ///< the size of the input alphabet
    long len_;            ///< the length of the sample sequence from which the estimates are calculated

    double nPseudoCounts_ = 1;
    double zConfidenceIntervalSize_ = 0;
    double minimumVariance_ = 0;
    double exponent_ = 1;

    SWIGCODE(%ignore State::operator=;)
    /**
     * Estimates for a sequence \f$x\f$ are computed by parsing the sequence from left to right, while gathering the required statistics (occurrence counts in the underlying sample sequence) by traversing the suffix tree accordingly.
     *
     * A `State` object captures the relevant information during the parsing of sequence \f$x\f$, and after finalizing the parsing by calling `eval()`.
     *
     * Therefore, typically do:
     *
     * > `state.reset();`
     * > `extendBy(sequence);`
     * > `extendBy(o, u);`
     * > `extendBy(sequence);`
     * > `eval();`
     *
     * The estimate is then stored in `state_.f_` and the corresponding variance estimate (if `estimateVariance_ = true`) is stored in `state_.v_`.
     */
    struct State {
        State(const std::shared_ptr<stree::STree> stree) : pos_(stree) { }

        State& operator=(const State& state) {
            pos_.set(state.pos_);
            f_ = state.f_;
            v_ = state.v_;
            fB_ = state.fB_;
            len_ = state.len_;
            return *this;
        }

        void reset() {
            pos_.setRoot();
            f_ = fB_ = v_ = 1;
            len_ = 0;
        }

        stree::Position pos_; ///< the position in the suffix tree for the currently estimated sequence.
        double f_ = 1;        ///< related to the current estimate (used internally in different ways)
        double v_ = 1;
        double fB_ = 1;       ///< related to the current Bayesian estimate (used internally in different ways)
        long len_ = 0;
    };

    State state_;                   ///< captures the relevant information during the parsing of a sequence
    bool estimateVariance_ = true;  ///< should a variance estimate also be computed?

    const std::shared_ptr<stree::STree> stree_; ///< the suffix tree for the underlying sample sequence

    void extendBy(Symbol o, Symbol u = 0) {
        state_.len_++;
        if (nU_ == 0) {
            double ch = state_.pos_.count();
            state_.pos_.toSymbol(o);
            double cz = state_.pos_.count();
            state_.fB_ *= (ch == 0 ? 0.5 : (cz + 0.5 * nPseudoCounts_) / (ch + nPseudoCounts_));
        } else { // IO-case
            state_.pos_.toSymbol(u);
            double cu = state_.pos_.count(); // (c)ount of input (u) after current history
            state_.pos_.toSymbol(o);
            double co = state_.pos_.count(); // (c)ount of output (o) after current history and input u
            if (cu == 0) state_.f_ *= (1.0 / nO_);
            else         state_.f_ *= (co / cu) ;

            if (estimateVariance_) {
                if (cu < 2) { state_.v_ *= 0.75; }
                else {
                    //double pB = (co + 1/nO_ * nPseudoCounts_) / (cu + nPseudoCounts_);
                    double pBv = (co + 0.5 * nPseudoCounts_) / (cu + nPseudoCounts_);
                    double vB = pBv * (1 - pBv) / (cu);
                    double confidenceIntervalRadius =
                            zConfidenceIntervalSize_ * std::sqrt(vB);
                    if (pBv <= 0.5) pBv = std::min(pBv + confidenceIntervalRadius, 0.5);
                    else pBv = std::max(0.5, pBv - confidenceIntervalRadius);
                    state_.fB_ *= pBv;
                    vB = pBv * (1 - pBv) / (cu);

                    state_.v_ *= std::max(pBv * pBv - vB, 0.0);
                }
            }
        }
    }

    void extendBy(const Sequence &seq) { for (long i = 0; i < seq.length(); ++i) { extendBy(seq.o(i), seq.u(i)); } }

    /** Call this after parsing a sequence by calls to `extendBy()` to finalize the parsing. This places the estimate in `state_.f_` and the corresponding variance estimate (if `estimateVariance = true`) in `state_.v_`.
     */
    void eval() {
        if (nU_ == 0) {
            double n = std::max(len_ - state_.len_ + 1, 1L);
            double x = state_.pos_.count();
            state_.f_ = x / n;

            if (estimateVariance_) {
                // Compute a confidence interval for the estimate f, and use the bound that is closer to 1/2 to estimate the variance:
                // double fB = (x + 0.5 * nPseudoCounts_) / (n + nPseudoCounts_);
                double fB = state_.fB_;
                double confidenceIntervalRadius =
                        zConfidenceIntervalSize_ * std::sqrt(fB * (1 - fB) / (n + nPseudoCounts_));
                if (fB <= 0.5) fB = std::min(fB + confidenceIntervalRadius, (double) (0.5));
                else fB = std::max((double) (0.5), fB - confidenceIntervalRadius);
                state_.v_ = fB * (1 - fB) / (n);
            }
        }
        else if (estimateVariance_) { // IO-case
            state_.v_ = std::max(state_.fB_ * state_.fB_ - state_.v_, 0.0);
        }
        if (estimateVariance_) {
            state_.v_ = std::pow(state_.v_, exponent_);
            if (state_.len_ == 0) {
                state_.v_ = 1.0 / (nO_ * (double)(len_));
            }
            if (minimumVariance_ == 1) {
                state_.v_ = std::max(state_.v_, 1.0/(nO_ * (double)len_ * (double)len_));
                state_.v_ = std::max(state_.v_, 1e-15);
            } else {
                state_.v_ = std::max(state_.v_, minimumVariance_);
            }
        }
    }

}; // class Estimator

} // namespace tom
