namespace tom {

/** Devide the given `matrix` by its element-sum, i.e., normalize the matrix to have an element-sum of one, and return the element-sum.
 */
template<typename T>
double normalize(const DenseBase<T> &matrix) {
    double mat_sum = matrix.sum();
    const_cast< DenseBase<T> & >(matrix) /= mat_sum;
    return mat_sum;
}
SWIGCODE(%template(normalize) normalize<MatrixMd>;)

/** Devide each column of the given `matrix` by its sum, i.e., normalize the columns to have column-sum one. Return `true` if successful, or `false` if a column could not be normalized due to a zero column-sum.
 */
template<typename T>
bool normalizeCols(const DenseBase<T> &matrix) {
    double col_sum;
    bool success = true;
    for (int j = 0; j < matrix.cols(); ++j) {
        col_sum = matrix.col(j).sum();
        if (col_sum == 0) { success = false; }
        const_cast< DenseBase<T> & >(matrix).col(j) /= col_sum;
    }
    return success;
}
SWIGCODE(%template(normalizeCols) normalizeCols<MatrixMd>;)

/** Devide each row of the given `matrix` by its sum, i.e., normalize the rows to have row-sum one. Return `true` if successful, or `false` if a row could not be normalized due to a zero row-sum.
 */
template<typename T>
bool normalizeRows(const DenseBase<T> &matrix) {
    double row_sum;
    bool success = true;
    for (int i = 0; i < matrix.rows(); ++i) {
        row_sum = matrix.row(i).sum();
        if (row_sum == 0) { success = false; }
        const_cast< DenseBase<T> & >(matrix).row(i) /= row_sum;
    }
    return success;
}
SWIGCODE(%template(normalizeRows) normalizeRows<MatrixMd>;)

/**
 * Return the Kronecker-product \f$A\otimes B\f$ of the matrices `A` and `B`. */
template< typename D1, typename D2 >
MatrixXd kron(const MatrixBase<D1>& A, const MatrixBase<D2>& B) {
    MatrixXd result(A.rows() * B.rows(), A.cols() * B.cols());
    for (long j = 0; j < A.cols(); ++j) {
        for (long i = 0; i < A.rows(); ++i) {
            result.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i,j) * B;
        }
    }
    return result;
}
SWIGCODE(%template(kron) kron<MatrixMd, MatrixMd>;)

SWIGCODE(%kwargs;)

/** Return the column-wise generalized mean with exponent `p` (default 1) of the given `matrix`. For `p` = 1, 0, -1 this is the arithmetic, geometric and harmonic mean, respectively.
 *
 * Note that for values of `p` other than {1, 2k} this requires all matrix entries to be positive.
 */
template<typename T>
RowVectorXd colwiseMean(const MatrixBase <T> &matrix, double p = 1.0) {
    RowVectorXd result(matrix.cols());
    if (p == std::numeric_limits<double>::infinity()) { result = matrix.array().colwise().maxCoeff(); }
    else if (p == - std::numeric_limits<double>::infinity()) { result = matrix.array().colwise().minCoeff(); }
    else if (p == 0) { // geometric mean
        result = matrix.array().abs().log().colwise().sum() / matrix.rows();
        result = result.array().exp();
    } else if (p == 1) { // arithmetic mean
        result = matrix.array().colwise().sum() / matrix.rows();
    } else {
        result = matrix.array().pow(p).colwise().sum() / matrix.rows();
        result = result.array().pow(1.0 / p);
    }
    return result;
}
SWIGCODE(%template(colwiseMean) colwiseMean<MatrixMd>;)

/** Return the row-wise generalized mean with exponent `p` (default 1) of the given `matrix`. For `p` = 1, 0, -1 this is the arithmetic, geometric and harmonic mean, respectively.
 *
 * Note that for values of `p` other than {1, 2k} this requires all matrix entries to be positive.
 */
template<typename T>
VectorXd rowwiseMean(const MatrixBase <T> &matrix, double p = 1.0) { return colwiseMean(matrix.transpose(), p).transpose(); }
SWIGCODE(%template(rowwiseMean) rowwiseMean<MatrixMd>;)

/** Return the weighted norm of `M` with weights given in `W`, or the squared weighted norm if `squared` is set to `true`.
 * Depending on the size of `W`, the given weights are interpreted in different ways, assuming `M` is of size m x n:
 *
 * - if `W` is of size zero, then no weights are used and the Frobenius norm |M|_F is computed
 * - if `W` is of size m+n x 1, then row and column weights [w_r; w_c] = W are assumed and |M|_D(w_r w_c^T) is computed
 * - if `W` is of size m x n, then element-wise weights are assumed and |M|_D(W) is computed
 * - if `W` is of size m x mn, then a block-diagonal weight matrix is assumed and |M|_D(W1,...,Wn) is computed
 * - if `W` is of size mn x mn, then a full weight matrix is assumed and |M|_W is computed
 */
template<typename T, typename T1>
double weightedNorm(const MatrixBase<T> &M, const MatrixBase<T1> &W, bool squared = false) throw(std::invalid_argument) {
    double result = 0;
    if (W.size() == 0) { result = M.squaredNorm(); }
    else if (W.cols() == 1 and W.rows() == M.rows() + M.cols()) {
        result = (W.col(0).head(M.rows()).asDiagonal() * M * W.col(0).tail(M.cols()).asDiagonal()).squaredNorm();
    } else if (W.rows() == M.rows()) {
        if (W.cols() == M.cols()) { result = (W.array() * M.array().square()).sum(); }
        else if (W.cols() == M.cols() * W.rows()) {
            for (long j = 0; j < M.cols(); ++j) {
                result += M.col(j).transpose() * W.middleCols(j * W.rows(), W.rows()) * M.col(j);
            }
        } else { throw std::invalid_argument("size mismatch for W and M"); }
    } else if (W.cols() == M.size() and W.rows() == M.size()) {
        MatrixXd Mcopy = M;
        Map<VectorXd, 1> vecM(Mcopy.data(), Mcopy.size());
        result = vecM.transpose() * W * vecM;
    } else { throw std::invalid_argument("size mismatch for W and M"); }
    result /= M.size();
    return squared ? result : std::sqrt(result);
}
SWIGCODE(%template(weightedNorm) weightedNorm<MatrixMd, MatrixMd>;)

SWIGCODE(%apply const MatrixBase<MatrixXd>& OUTPUT { const MatrixBase<MatrixXd>& X };)

//<editor-fold desc="Solve implementations">
/**
 * Return
 * \ifnot PY
 * in the output-argument `X`
 * \endif
 * the ordinary least-squares (OLS) solution to the problem `A` * `X` = `M` (or if `transposed` to `X` * `A` = `M`) using a `method` from {"Cholesky", "LDLT", "QR" (default), "SVD", "JacobiSVD"}.
 *
 * The "Cholesky" method solves the normal equations using a Cholesky decomposition. This is the fastest method, but loses most precision and requires the problem to be overdetermined and `A` to have full rank.
 *
 * The "LDLT" method is essentially the same as "Cholesky", but uses a more robust Cholesky decomposition with pivoting that also avoids taking a square root. This method is recommended over "Cholesky" by Eigen3.
 *
 * The "QR" method uses a QR decomposition. This is slower than "Cholesky", but gives more precision. The marix `A` should have full rank.
 *
 * The "SVD" uses an SVD decomposition. This is the slowest, but gives best precision. Also, the matrix `A` does not need to have full rank, and in the case of an underdetermined problem, the least-squares solution with the smallest norm is returned.
 *
 * The "JacobiSVD" method is similar to the "SVD" method, but uses a different (slower, but potentially more accurate) svd algorithm.
 */
template< typename D, typename D1, typename D2 >
C1(void) PY1(MatrixXd)
solveOLS(C2(const MatrixBase<D> &X,) const MatrixBase<D1> &A, const MatrixBase<D2> &M, bool transposed = false, const std::string& method = "QR") {
    if (transposed) {
        if (method == "Cholesky") {
            const_cast<MatrixBase<D> &>(X).derived().noalias() = (A * A.transpose()).llt().solve(A * M.transpose()).transpose();
        } else if (method == "LDLT") {
            const_cast<MatrixBase<D> &>(X).derived().noalias() = (A * A.transpose()).ldlt().solve(A * M.transpose()).transpose();
        } else if (method == "QR") {
            const_cast<MatrixBase<D> &>(X).derived().noalias() = A.transpose().colPivHouseholderQr().solve(M.transpose()).transpose();
        } else if (method == "SVD") {
            const_cast<MatrixBase<D> &>(X).derived().noalias() = A.transpose().bdcSvd(ComputeThinU | ComputeThinV).solve(M.transpose()).transpose();
        } else if (method == "JacobiSVD") {
            const_cast<MatrixBase<D> &>(X).derived().noalias() = A.transpose().jacobiSvd(ComputeThinU | ComputeThinV).solve(M.transpose()).transpose();
        } else { throw std::invalid_argument("unrecognized method"); }
    } else {
        if (method == "Cholesky") {
            const_cast<MatrixBase<D> &>(X).derived().noalias() = (A.transpose() * A).llt().solve(A.transpose() * M);
        } else if (method == "LDLT") {
            const_cast<MatrixBase<D> &>(X).derived().noalias() = (A.transpose() * A).ldlt().solve(A.transpose() * M);
        } else if (method == "QR") {
            const_cast<MatrixBase<D> &>(X).derived().noalias() = A.colPivHouseholderQr().solve(M);
        } else if (method == "SVD") {
            const_cast<MatrixBase<D> &>(X).derived().noalias() = A.bdcSvd(ComputeThinU | ComputeThinV).solve(M);
        } else if (method == "JacobiSVD") {
            const_cast<MatrixBase<D> &>(X).derived().noalias() = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(M);
        } else { throw std::invalid_argument("unrecognized method"); }
    }
}
SWIGCODE(%template(solveOLS) solveOLS<MatrixXd, MatrixMd, MatrixMd>;)

/**
 * Return
 * \ifnot PY
 * in the output-argument `X`,
 * \endif
 * the row or column weighted least-squares solution to the problem `A` * `X` = `M` with row-weights given in the column vector `W` (or if `transposed` to `X` * `A` = `M` with column-weights given in the row vector `W`) using a `method` from {"Cholesky", "LDLT" (default), "QR", "SVD", "JacobiSVD"}.
 *
 * This computes `X` that minimizes |D(sqrt_W) * (`A` * `X` - `M`)|_F (or |(`X` * `A` - `M`) * D(sqrt_W)|_F if `transposed`), where `sqrt_W` is the element-wise square-root of `W`, i.e., `W` = `sqrt_W` .* `sqrt_W`, and `.*` denotes the element-wise product. The computation is done by reducing the problem to an OLS problem that is then solved according to the given `method` as detailed below (see also `solveOLS()`). Note that the weights in `W` must be strictly greater than zero.
 *
 * Note that column weights have no effect in the default case, and row weight have no effect if `transposed`, and are therefore ommitted.
 *
 * The "Cholesky" method solves the normal equations using a Cholesky decomposition. This is the fastest method, but loses most precision and requires the problem to be overdetermined and `A` to have full rank.
 *
 * The "LDLT" method is essentially the same as "Cholesky", but uses a more robust Cholesky decomposition with pivoting that also avoids taking a square root. This method is recommended over "Cholesky" by Eigen3.
 *
 * The "QR" method uses a QR decomposition. This is slower than "Cholesky", but gives more precision. The marix `A` should have full rank.
 *
 * The "SVD" uses an SVD decomposition. This is the slowest, but gives best precision. Also, the matrix `A` does not need to have full rank, and in the case of an underdetermined problem, the least-squares solution with the smallest norm is returned.
 *
 * The "JacobiSVD" method is similar to the "SVD" method, but uses a different (slower, but potentially more accurate) svd algorithm.
 */
template< typename D, typename D1, typename D2, typename D3>
C1(void) PY1(MatrixXd)
solveRowColWLS(C2(const MatrixBase<D> &X,) const MatrixBase<D1>&A, const MatrixBase<D2>& M, const MatrixBase<D3>& W, bool transposed = false, const std::string &method = "LDLT") {
    if (transposed) {
        if (method == "Cholesky") {
            const_cast<MatrixBase<D> &>(X).derived().noalias() = (A * W.asDiagonal() * A.transpose()).llt().solve(
                    A * W.asDiagonal() * M.transpose()).transpose();
        } else if (method == "LDLT") {
            const_cast<MatrixBase<D> &>(X).derived().noalias() = (A * W.asDiagonal() * A.transpose()).ldlt().solve(
                    A * W.asDiagonal() * M.transpose()).transpose();
        } else {
            RowVectorXd sqrt_W = W.cwiseSqrt();
            solveOLS(X, A * sqrt_W.asDiagonal(), M * sqrt_W.asDiagonal(), transposed, method);
        }
    } else {
        if (method == "Cholesky") {
            const_cast<MatrixBase<D> &>(X).derived().noalias() = (A.transpose() * W.asDiagonal() * A).llt().solve(
                    A.transpose() * W.asDiagonal() * M);
        } else if (method == "LDLT") {
            const_cast<MatrixBase<D> &>(X).derived().noalias() = (A.transpose() * W.asDiagonal() * A).ldlt().solve(
                    A.transpose() * W.asDiagonal() * M);
        } else {
            VectorXd sqrt_W = W.cwiseSqrt();
            solveOLS(X, sqrt_W.asDiagonal() * A, sqrt_W.asDiagonal() * M, transposed, method);
        }
    }
}
SWIGCODE(%template(solveRowColWLS) solveRowColWLS<MatrixXd, MatrixMd, MatrixMd, MatrixMd>;)

/**
 * Return
 * \ifnot PY
 * in the output-argument `X`,
 * \endif
 * the (element-wise) D(`W`)-weighted least-squares (WLS) solution to the problem `A` * `X` = `M` (or to `X` * `A` = `M` if `transposed`) using a `method` from {"Cholesky", "LDLT" (default), "QR", "SVD", "JacobiSVD"}.
 *
 * This computes `X` that minimizes |`A` * `X` - `M`|_D(`W`) (or |`X` * `A` - `M`|_D(`W`) if `transposed`). The computation is done by reducing the problem to a set of OLS problems that are then solved according to the given `method` as detailed below (see also `solveOLS()`). Note that the weights in `W` must be strictly greater than zero.
 *
 * The "Cholesky" method solves the normal equations using a Cholesky decomposition. This is the fastest method, but loses most precision and requires the problem to be overdetermined and `A` to have full rank.
 *
 * The "LDLT" method is essentially the same as "Cholesky", but uses a more robust Cholesky decomposition with pivoting that also avoids taking a square root. This method is recommended over "Cholesky" by Eigen3.
 *
 * The "QR" method uses a QR decomposition. This is slower than "Cholesky", but gives more precision. The marix `A` should have full rank.
 *
 * The "SVD" uses an SVD decomposition. This is the slowest, but gives best precision. Also, the matrix `A` does not need to have full rank, and in the case of an underdetermined problem, the least-squares solution with the smallest norm is returned.
 *
 * The "JacobiSVD" method is similar to the "SVD" method, but uses a different (slower, but potentially more accurate) svd algorithm.
 */
template< typename D, typename D1, typename D2, typename D3>
C1(void) PY1(MatrixXd)
solveWLS(C2(const MatrixBase<D> &X,) const MatrixBase<D1>& A, const MatrixBase<D2>& M, const MatrixBase<D3>& W, bool transposed = false, const std::string &method = "LDLT") {
    if (transposed) {
        const_cast<MatrixBase<D> &>(X).derived().resize(M.rows(), A.rows());
        if (method == "Cholesky") {
            LLT<MatrixXd> llt;
            MatrixXd A_WMT = A * (W.array() * M.array()).matrix().transpose();
            for (int i = 0; i < M.rows(); ++i) {
                llt.compute(A * W.row(i).asDiagonal() * A.transpose());
                const_cast<MatrixBase<D> &>(X).row(i).noalias() = llt.solve(A_WMT.col(i)).transpose();
            }
        } else if (method == "LDLT") {
            LDLT<MatrixXd> ldlt;
            MatrixXd A_WMT = A * (W.array() * M.array()).matrix().transpose();
            for (int i = 0; i < M.rows(); ++i) {
                ldlt.compute(A * W.row(i).asDiagonal() * A.transpose());
                const_cast<MatrixBase<D> &>(X).row(i).noalias() = ldlt.solve(A_WMT.col(i)).transpose();
            }
        } else {
            MatrixXd sqrt_W = W.cwiseSqrt();
            for (int i = 0; i < M.rows(); ++i) {
                solveOLS(const_cast<MatrixBase<D> &>(X).row(i), A * sqrt_W.row(i).asDiagonal(), M.row(i) * sqrt_W.row(i).asDiagonal(), transposed, method);
            }
        }
    } else {
        const_cast<MatrixBase<D> &>(X).derived().resize(A.cols(), M.cols());
        if (method == "Cholesky") {
            LLT<MatrixXd> llt;
            MatrixXd AT_WM = A.transpose() * (W.array() * M.array()).matrix();
            for (int j = 0; j < M.cols(); ++j) {
                llt.compute(A.transpose() * W.col(j).asDiagonal() * A);
                const_cast<MatrixBase<D> &>(X).col(j).noalias() = llt.solve(AT_WM.col(j));
            }
        } else if (method == "LDLT") {
            LDLT<MatrixXd> ldlt;
            MatrixXd AT_WM = A.transpose() * (W.array() * M.array()).matrix();
            for (int j = 0; j < M.cols(); ++j) {
                ldlt.compute(A.transpose() * W.col(j).asDiagonal() * A);
                const_cast<MatrixBase<D> &>(X).col(j).noalias() = ldlt.solve(AT_WM.col(j));
            }
        } else {
            MatrixXd sqrt_W = W.cwiseSqrt();
            for (int j = 0; j < M.cols(); ++j) {
                solveOLS(const_cast<MatrixBase<D> &>(X).col(j), sqrt_W.col(j).asDiagonal() * A, sqrt_W.col(j).asDiagonal() * M.col(j), transposed, method);
            }
        }
    }
}
SWIGCODE(%template(solveWLS) solveWLS<MatrixXd, MatrixXd, MatrixXd, MatrixXd>;)

/**
 * Return
 * \ifnot PY
 * in the output-argument `X`,
 * \endif
 * the D(W1,..., Wm)-weighted least-squares (GLS) solution to the overdetermined problem `A` * `X` = `M` (or to `X` * `A` = `M` if `transposed`) using a `method` from {"Cholesky", "LDLT" (default)}, where the block-diagonal symmetric and positive definite weight matrix is given by `W` = [W1,..., Wn], where each `Wj` is the full weight matrix for the column j of `M`.
 *
 * This computes `X` that minimizes |`A` * `X` - `M`|_D(W1,...,Wn) (or |`X` * `A` - `M`|_D(W1,...,Wn) if `transposed`).
 *
 * Note that the "LDLT" method is essentially the same as "Cholesky", but uses a more robust Cholesky decomposition with pivoting that also avoids taking a square root. This method is recommended over "Cholesky" by Eigen3.
 */
template< typename D, typename D1, typename D2, typename D3>
C1(void) PY1(MatrixXd)
solveGLS(C2(const MatrixBase<D> &X,) const MatrixBase<D1> &A, const MatrixBase<D2>& M, const MatrixBase<D3>& W, bool transposed = false, const std::string &method = "LDLT") {
    assert(W.rows() == M.rows() and W.cols() == M.cols() * M.rows()); // Block-diagonal weight matrix
    if (transposed) {
        // solve XA≈M by solving [\sum_j AjAj^T ⊗ Wj] vec(X) = vec( [W1M1, W2M2, ..., WnMn] A^T )
        MatrixXd AtI_W_ATtI(A.rows() * W.rows(), A.rows() * W.rows()); // we only compute lower triagonal part
        MatrixXd AiijAjjjWj(W.rows(), W.rows());
        for (long jj = 0; jj < A.rows(); ++jj) {
            for (long ii = jj; ii < A.rows(); ++ii) {
                AiijAjjjWj.setZero(W.rows(), W.rows());
                for (long j = 0; j < M.cols(); ++j) {
                    AiijAjjjWj += A(ii, j) * A(jj, j) * W.middleCols(j * W.rows(), W.rows());
                }
                AtI_W_ATtI.block(ii * W.rows(), jj * W.rows(), W.rows(), W.rows()) = AiijAjjjWj;
            }
        }
        MatrixXd WM = MatrixXd(W.rows(), M.cols());
        for (int j = 0; j < M.cols(); ++j) {
            WM.col(j).noalias() = W.middleCols(j*W.rows(), W.rows()) * M.col(j);
        }
        WM *= A.transpose();
        const Map<const VectorXd> vecWMAT(WM.data(),WM.size());
        MatrixXd vecX;
        if (method == "LDLT") { vecX.noalias() = AtI_W_ATtI.ldlt().solve(vecWMAT); }
        else if (method == "Cholesky") { vecX.noalias() = AtI_W_ATtI.llt().solve(vecWMAT); }
        else { throw std::invalid_argument("unrecognized method"); }
        vecX.resize(M.rows(), A.rows());
        const_cast<MatrixBase<D>&>(X).derived() = std::move(vecX);
    } else {
        // solve AX≈M by solving SjAXj ≈ SjMj for each column j
        if (method == "Cholesky") {
            LLT<MatrixXd> llt;
            const_cast<MatrixBase<D> &>(X).derived().resize(A.cols(), M.cols());
            for (int j = 0; j < X.cols(); ++j) {
                llt.compute(A.transpose() * W.middleCols(j * W.rows(), W.rows()) * A);
                const_cast<MatrixBase<D> &>(X).col(j).noalias() = llt.solve(
                        A.transpose() * (W.middleCols(j * W.rows(), W.rows()) * M.col(j)));
            }
        } else if (method == "LDLT") {
            LDLT<MatrixXd> ldlt;
            const_cast<MatrixBase<D> &>(X).derived().resize(A.cols(), M.cols());
            for (int j = 0; j < X.cols(); ++j) {
                ldlt.compute(A.transpose() * W.middleCols(j * W.rows(), W.rows()) * A);
                const_cast<MatrixBase<D> &>(X).col(j).noalias() = ldlt.solve(
                        A.transpose() * (W.middleCols(j * W.rows(), W.rows()) * M.col(j)));
            }
        }
    }
}
SWIGCODE(%template(solveGLS) solveGLS<MatrixXd, MatrixXd, MatrixXd, MatrixXd>;)
//</editor-fold>

/**
 * Return
 * \ifnot PY
 * in the output-argument `X`,
 * \endif
 * the generalized weighted least-squares (GLS) solution to the overdetermined problem `A` * `X` = `M` (or to the problem `X` * `A` = `M` if `transposed`) with the given weights `W` using a `method` from {"Cholesky", "LDLT" (default), "QR", "SVD", "JacobiSVD"}.
 *
 * This computes `X` that minimizes the weighted norm |`A` * `X` - `M`|_W (or |`X` * `A` - `M`|_W if `transposed`), utilizing the structure of W that depends on the size of the supplied weights `W` in the following way, assuming `M` is of size m x n:
 *
 * - If `W` is of size zero (default), then no weights are used and the ordinary least squares (OLS) solution minimizing the Frobenius norm is returned.
 * - If `W` is of size m x n, then element-wise weights are assumed, i.e., W = D(`W`), resulting in the weighted least squares (WLS) solution.
 * - If `W` is of size m x nm, then a block-diagonal structure for W is assumed, i.e., W = D(`W`_1, ..., `W`_m), where `W`_j is the j-th (m x m)-block of `W` corresponding to the weights for [`M`]_j, which must be symmetric an positive definite.
 * - If `W` is of size m x 1, then these are regarded as row-weights, which only make sense if not `transposed`.
 * - If `W` is of size 1 x n, then these are regarded as column-weights, which only make sense if `transposed`.
 * - If `W` is of size m+n x 1 then `W`[:m] are regarded as row-weights and `W`[m:] as column-weights. If `transposed` the row-weights, else the column-weights, have no effect.
 *
 * Note that solving GLS with full weight matrix is expensive, and therefore only solving with block-diagonal structured weight matrix is supported, and then only the methods "Cholesky" and "LDLT" which can make use of extra simplifications.
 *
 * The computation is done by reducing the problem to a set of OLS problems that are then solved according to the given `method` as detailed below. Note that the weights in `W` must be strictly greater than zero.
 *
 * The "Cholesky" method solves the normal equations using a Cholesky decomposition. This is the fastest method, but loses most precision and requires the problem to be overdetermined and `A` to have full rank.
 *
 * The "LDLT" method is essentially the same as "Cholesky", but uses a more robust Cholesky decomposition with pivoting that also avoids taking a square root. This method is recommended over "Cholesky" by Eigen3.
 *
 * The "QR" method uses a QR decomposition. This is slower than "Cholesky", but gives more precision. The marix `A` should have full rank.
 *
 * The "SVD" uses an SVD decomposition. This is the slowest, but gives best precision. Also, the matrix `A` does not need to have full rank in which case the least-squares solution with the smallest norm is returned.
 *
 * The "JacobiSVD" method is similar to the "SVD" method, but uses a different (slower, but potentially more accurate) svd algorithm.
 *
 */
template< typename D>
C1(void) PY1(MatrixXd)
solveLS(C2(const MatrixBase<D> &X,) const MatrixXd &A, const MatrixXd &M, const MatrixXd &W = MatrixXd(), bool transposed = false, std::string method = "LDLT") throw(std::invalid_argument) {
    if (not (method == "Cholesky" or method == "LDLT" or method == "QR" or method == "SVD" or method == "JacobiSVD")) throw std::invalid_argument("method not recognized");
    if (not transposed) { // solve AX=M
        if (A.rows() != M.rows()) throw std::invalid_argument("LHS rows != RHS rows");
        if (A.cols() > A.rows()) throw std::invalid_argument("system must be overdetermined");
        if (W.size() == 0) { solveOLS(X, A, M, transposed, method); }
        else if (W.rows() == M.rows() and W.cols() == M.cols()) { solveWLS(X, A, M, W, transposed, method); }
        else if (W.rows() == M.rows() and W.cols() == M.rows() * M.cols()) {
            if (method != "Cholesky" and method != "LDLT") throw std::invalid_argument("GLS only supported with method 'Cholesky' or 'LDLT'");
            solveGLS(X, A, M, W, transposed, method); }
        else if (W.rows() == W.cols() and W.cols() == M.rows() * M.cols()) throw std::invalid_argument("GLS with full weight matrix not implemented");
        else if (W.rows() == M.rows() and W.cols() == 1) { solveRowColWLS(X, A, M, W, transposed, method); }
        else if (W.cols() == 1 and W.rows() == M.rows() + M.cols()) { solveRowColWLS(X, A, M, W.col(0).head(M.rows()), transposed, method); }
        else { throw std::invalid_argument("size mismatch for M and W"); }
    } else { // solve XA=M
        if (A.cols() != M.cols()) throw std::invalid_argument("LHS cols != RHS cols");
        if (A.rows() > A.cols()) throw std::invalid_argument("system must be overdetermined");
        if (W.size() == 0) { solveOLS(X, A, M, transposed, method); }
        else if (W.rows() == M.rows() and W.cols() == M.cols()) { solveWLS(X, A, M, W, transposed, method); }
        else if (W.rows() == M.rows() and W.cols() == M.rows() * M.cols()) {
            if (method != "Cholesky" and method != "LDLT") throw std::invalid_argument("GLS only supported with method 'Cholesky' or 'LDLT'");
            solveGLS(X, A, M, W, transposed, method); }
        else if (W.rows() == W.cols() and W.cols() == M.rows() * M.cols()) throw std::invalid_argument("LS with full weight matrix not implemented");
        else if (W.rows() == 1 and W.cols() == M.cols()) { solveRowColWLS(X, A, M, W, transposed, method); }
        else if (W.cols() == 1 and W.rows() == M.rows() + M.cols()) { solveRowColWLS(X, A, M, W.col(0).tail(M.cols()).transpose(), transposed, method); }
        else { throw std::invalid_argument("size mismatch for M and W"); }
    }
}
SWIGCODE(%template(solveLS) solveLS<MatrixXd>;)

/**
 * Return a new weight matrix for `X`, assuming `X` is a solution to the D(`W`)-weighted WLS problem `B` * `X` = `M`.
 *
 * Note that the columns of `X` can be regarded as coordinate representations for the columns of `M` with respect to a basis given by the columns of `B`. This function transforms the given weights for the columns of `M` to appropriate weights for the coordinates in the columns of `X`. The resulting weight matrix for `X` will therefore be block-diagonal in general, but if `covariances` is set to `false`, the off-diagonal weights are ignored, resulting in element-wise weights for `X`.
 *
 * The returned matrix will therefore be
 *
 * - [B^T * D([W]_1) * B, ..., B^T * D([W]_m) * B] of size B.cols() x B.cols * M.cols() if `covariances`
 * - [diag(B^T * D([W]_1) * B), ..., diag(B^T * D([W]_m) * B)] of size B.cols() x M.cols() otherwise.
 */
template< typename D1, typename D2 >
MatrixXd transformWeights(const MatrixBase<D1> &W, const MatrixBase<D2> &B, bool covariances = true) {
    MatrixXd newWeights(B.cols(), W.cols() * (covariances ? B.cols() : 1));
    for (int j = 0; j < W.cols(); ++j) {
        if (covariances) {
            newWeights.middleCols(j * B.cols(), B.cols()).noalias() = B.transpose() * W.col(j).asDiagonal() * B;
        } else {
            newWeights.col(j).noalias() = (B.transpose() * W.col(j).asDiagonal() * B).diagonal();
        }
    }
    return newWeights;
}
SWIGCODE(%template(transformWeights) transformWeights<MatrixMd, MatrixMd>;)

/**
 * Return the pseudo-inverse of the given matrix `M` computed according to the given `method` from {"Cholesky", "QR", "SVD" (default), "JacobiSVD"}.
 *
 * If `method` is "SVD" or "JacobiSVD", the classical pseudo-inverse is computed from the svd of `M`.
 *
 * If `method` is "QR", the pseudo-inverse is computed from the QR-factorization of `M`, which requires `M` to have full rank.
 *
 * If `method` is "Cholesky" or "LDLT", the pseudo-inverse is computed as \f$(M^\top M)^{-1} M^\top\f$ or \f$M^\top (M M^\top)^{-1}\f$ depending on the size of `M`, which requires `M` to have full rank.
 */
template<typename D>
MatrixXd pinv(const MatrixBase<D> &M, const std::string &method = "SVD") {
    MatrixXd result;
    unsigned long I_size = M.rows();
    bool transpose = false;
    if ((method == "LDLT" or method == "Cholesky" or method == "QR") and (M.rows() < M.cols())) {
        I_size = M.cols();
        transpose = true;
    }
    solveOLS(result, M, MatrixXd::Identity(I_size, I_size), transpose, method);
    return result;
}
SWIGCODE(%template(pinv) pinv<MatrixMd>;)

/**
 * Compute in the arguments `B` and `A` (an approximation to) the best weighted rank-d approximation to `M` with element-wise weights `W` such that |`B` * `A` - `M`|_D(W) is minimized. This is computed iteratively starting from an initial approximation given by `B` * `A` using "alternating projections" solved via the given `method`. The termination of the iteration is controlled by the given `stopCondition`.
 * \ifnot PY
 * See `StopCondition`.
 * \else
 * See `tom.util.StopCondition`.
 * \endif
 */
template<typename D1, typename D2, typename D3, typename D4>
void improveWLRA(const MatrixBase<D1> &B, const MatrixBase<D2> &A, const MatrixBase<D3> &M, const MatrixBase<D4> &W,
            const StopCondition& stopCondition = StopCondition(50, 1e-5, 1e-12), const std::string &method = "Cholesky") {
    const_cast<StopCondition&>(stopCondition).reset();
    while (not const_cast<StopCondition&>(stopCondition)(weightedNorm(M - B * A, W))) {
        solveWLS(B, A, M, W, true, method);
        solveWLS(A, B, M, W, false, method);
    }
}
SWIGCODE(%template(improveWLRA) improveWLRA<MatrixMd, MatrixMd, MatrixXd, MatrixXd >;)

/**
 * \ifnot PY
 * Compute in the arguments `B` and `A`
 * \else
 * Return in a tuple [B,A]
 * \endif
 * (an approximation to) the best weighted rank-d approximation to `M` with element-wise weights `W` such that |`B` * `A` - `M`|_D(vect(W)) is minimized. This is computed iteratively starting from an initial approximation given by `B_init` using "alternating projections", which are in turn solved via the given `method`. The termination of the iteration is controlled by the given `stopCondition`.
 * \ifnot PY
 * See `StopCondition`.
 * \else
 * See `tom.util.StopCondition`.
 * \endif
 */
SWIGCODE(%apply const MatrixBase<MatrixXd>& OUTPUT { const MatrixBase<MatrixXd> &B, const MatrixBase<MatrixXd> &A };)
template<typename D1, typename D2, typename D3, typename D4, typename D5>
C1(void) PY2(tuple<MatrixXd, MatrixXd>)
computeWLRA(C3(const MatrixBase<D1> &B, const MatrixBase<D2> &A,) const MatrixBase<D3> &M, const MatrixBase<D4> &W, const MatrixBase<D5> &B_init, const StopCondition& stopCondition = StopCondition(50, 1e-5, 1e-12), const std::string &method = "Cholesky") throw(std::invalid_argument) {
    if (W.rows() != M.rows() or W.cols() != M.cols() or B_init.rows() != M.rows() or B_init.cols() > M.cols()) {
        throw std::invalid_argument("size mismatch");
    }
    const_cast<MatrixBase<D1>&>(B).derived() = B_init;
    solveWLS(A, B, M, W, false, method);
    improveWLRA(B, A, M, W, stopCondition, method);
}
SWIGCODE(%template(computeWLRA) computeWLRA<MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd >;)

SWIGCODE(%clear const MatrixBase<MatrixXd> &X, const MatrixBase<MatrixXd> &B, const MatrixBase<MatrixXd> &A;)
SWIGCODE(%clearkwargs;)

} // namespace tom
