namespace tom {

SWIGCODE(%feature("python:slot", "tp_repr", functype="reprfunc") Estimator::repr;)

/**
 * This class computes estimates for \f$f( x )\f$ and corresponding variance estimates for sequences \f$x\f$ based on a suffix tree representation of a sample sequence.
 */
class Estimator {
public:
    /**
     * Create an `Estimator` for a sample sequence data given by a suffix tree representation `stree`. */
    Estimator(const std::shared_ptr<stree::STree> &stree) : state_(stree), stree_(stree) {
        nO_ = stree_->sequence().nOutputSymbols();
        nU_ = stree_->sequence().nInputSymbols();
        N_ = (stree_->sequence().length());
        state_.pos_.setRoot();
        regularization("default");
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
        state_.reset();
        extendBy(sequence);
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
        state_.reset();
        extendBy(sequence);
        return this->v();
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
        state_.reset();
        extendBy(sequence);
        f = state_.f_;
        v = this->v();
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
                F.coeffRef(i, j) = state_.f_;
                V.coeffRef(i, j) = this->v();
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
    SWIGCODE(%apply double *OUTPUT { double *zConfidenceIntervalSizeOut, double *minimumVarianceOut, double *exponentOut };)
    SWIGCODE(%apply std::string *OUTPUT { std::string *confidenceIntervalTypeOut };)
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
    regularization(C5(double *zConfidenceIntervalSizeOut, std::string *confidenceIntervalTypeOut, double *minimumVarianceOut, double *exponentOut,)
                   double zConfidenceIntervalSize = -1,
                   std::string confidenceIntervalType = "",
                   double minimumVariance = -1,
                   double exponent = -1,
                   std::string preset = "") throw(std::invalid_argument) {
        if (preset == "default") {
            zConfidenceIntervalSize_ = 1;
            confidenceIntervalType_ = "Wilson";
            minimumVariance_         = 0;
            exponent_                = 1;
        } else if (preset == "none") {
            zConfidenceIntervalSize_ = 0;
            confidenceIntervalType_ = "Wilson";
            minimumVariance_         = 0;
            exponent_                = 1;
        } else if (preset != "") {
            throw std::invalid_argument("unrecognized preset");
        }
        if (zConfidenceIntervalSize != -1) this->zConfidenceIntervalSize_ = zConfidenceIntervalSize;
        if (confidenceIntervalType != "") this->confidenceIntervalType_ = confidenceIntervalType;
        if (minimumVariance         != -1) this->minimumVariance_ = minimumVariance;
        if (exponent                != -1) this->exponent_ = exponent;
        if (zConfidenceIntervalSizeOut) *zConfidenceIntervalSizeOut = this->zConfidenceIntervalSize_;
        if (confidenceIntervalTypeOut) *confidenceIntervalTypeOut = this->confidenceIntervalType_;
        if (minimumVarianceOut) *minimumVarianceOut = this->minimumVariance_;
        if (exponentOut) *exponentOut = this->exponent_;
    }

    SWIGCODE(%clear double *nPseudoCountsOut, double *zConfidenceIntervalSizeOut, double *minimumVarianceOut, double *exponentOut;)
    SWIGCODE(%clear std::string *confidenceIntervalTypeOut;)

#ifndef SWIG
#ifndef PY
    /**
     * Set the regularization parameters for the `Estimator`. This is a convenience function that simply calls `regularization(NULL, NULL, NULL, NULL, nPseudoCounts, zConfidenceIntervalSize, minimumVariance, exponent, preset);` */
    void regularization(double zConfidenceIntervalSize = -1,
                        std::string confidenceIntervalType = "",
                        double minimumVariance = -1,
                        double exponent = -1,
                        std::string preset = "") throw(std::invalid_argument) {
        regularization(NULL, NULL, NULL, NULL, zConfidenceIntervalSize, confidenceIntervalType, minimumVariance, exponent, preset);
    }
    /**
    * Set the regularization parameters for the `Estimator` to the given `preset`. This is a convenience function that simply calls `regularization(NULL, NULL, NULL, NULL, -1, -1, -1, -1, preset);` */
    void regularization(std::string preset) throw(std::invalid_argument) {
        regularization(NULL, NULL, NULL, NULL, -1, "", -1, -1, preset);
    }
#endif
#endif

    std::string confidenceIntervalType() const { return confidenceIntervalType_; }

    double zConfidenceIntervalSize() const { return zConfidenceIntervalSize_; }

    double minimumVariance() const { return minimumVariance_; }

    double exponent() const { return exponent_; }
    //</editor-fold>

private:
    int nO_;              ///< the size of the output alphabet
    int nU_;              ///< the size of the input alphabet
    long N_;            ///< the length of the sample sequence from which the estimates are calculated

    std::string confidenceIntervalType_;
    double zConfidenceIntervalSize_;
    double minimumVariance_;
    double exponent_;

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
            f2_ = state.f2_;
            v_ = state.v_;
            k_ = state.k_;
            return *this;
        }

        void reset() {
            pos_.setRoot();
            f_ = f2_ = 1;
            v_ = 0;
            k_ = 0;
        }

        stree::Position pos_; ///< the position in the suffix tree for the currently estimated sequence.
        double f_ = 1;        ///< the current estimate
        double f2_ = 1;
        double v_ = 0;
        long k_ = 0;          ///< the length of the currently parsed sequence
    };

    State state_;                   ///< captures the relevant information during the parsing of a sequence
    bool estimateVariance_ = true;  ///< should a variance estimate also be computed?

    const std::shared_ptr<stree::STree> stree_; ///< the suffix tree for the underlying sample sequence

    void extendBy(Symbol o, Symbol u = 0) {
        state_.k_++;
        /* `n` is the number of possible occurrences for `x_k`, and `c` is the number of occurrences of `x_k`. */
        if (nU_ != 0) { state_.pos_.toSymbol(u); }
        double n = state_.pos_.count();
        if (state_.pos_.isSuffix()) { --n; }
        state_.pos_.toSymbol(o);
        double c = state_.pos_.count();
        double p = (n == 0 ? 1.0 / nO_ : c / n);
        state_.f_ *= p;

        /* Next we compute a confidence interval for p: */
        double p_mid = (n == 0 ? 0.5 : c / n);
        double ci_lower = 0.0, ci_upper = 0.0;
        if (zConfidenceIntervalSize_ != 0) {
            double z2 = zConfidenceIntervalSize_ * zConfidenceIntervalSize_;
            ci_lower = ci_upper = 1.0;
            if (n > 0) {
                p_mid = (c + 0.5 * z2) / (n + z2);
                if (confidenceIntervalType_ == "Agresti–Coull") {
                    ci_lower = ci_upper = zConfidenceIntervalSize_ * std::sqrt(p_mid * (1 - p_mid) / (n + z2));
                } else if (confidenceIntervalType_ == "Wilson") {
                    ci_lower = ci_upper = zConfidenceIntervalSize_ * std::sqrt(n * p * (1 - p) + z2 / 4) / (n + z2);
                } else if (confidenceIntervalType_ == "Wilson_CC") { // Wilson score interval with continuity correction (more conservative):
                    if (c > 0) ci_lower = (0.5 + zConfidenceIntervalSize_ * std::sqrt(n * p * (1 - p) + z2 / 4 - 0.25 / n + p - 0.5)) / (n + z2);
                    if (c < n) ci_upper = (0.5 + zConfidenceIntervalSize_ * std::sqrt(n * p * (1 - p) + z2 / 4 - 0.25 / n - p + 0.5)) / (n + z2);
                }
            }
        }
        double p_lower = std::max(p_mid - ci_lower, 0.0);
        double p_upper = std::min(p_mid + ci_upper, 1.0);
        double p_centered = (p_mid < 0.5 ? std::min(p_upper, 0.5) : std::max(p_lower, 0.5));
        double p_extreme = std::min(p_lower, 1 - p_upper);

        double E_p2_overestimate = p_upper * p_upper;
        double var_p_overestimate = 0.5;
        double var_p_underestimate = 0.0;
        if (n > 0) {
            var_p_overestimate = p_centered * (1 - p_centered) / n;
            var_p_underestimate = p_extreme * (1 - p_extreme) / n;
        }

        state_.v_ = state_.v_ * (E_p2_overestimate - var_p_underestimate) + state_.f2_ * var_p_overestimate;
        state_.f2_ *= E_p2_overestimate;
    }

    void extendBy(const Sequence &seq) { for (long i = 0; i < seq.length(); ++i) { extendBy(seq.o(i), seq.u(i)); } }

    double v() const {
        //double var = state_.fB_ * (state_.fB_ - state_.fBm1_);
        double var = state_.v_;
        if (state_.k_ == 0) { var = 1e-15; }
        var = std::pow(var, exponent_);
        return var;
    }

}; // class Estimator

} // namespace tom
