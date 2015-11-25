namespace tom {

template<typename T>
double normalize(const DenseBase<T> &mat) {
    double mat_sum = mat.sum();
    if (mat_sum == 0) { return 0; }
    const_cast< DenseBase<T> & >(mat) /= mat_sum;
    return mat_sum;
}
SWIGCODE(%template(normalize) normalize<MatrixMd>;)

template<typename T>
bool normalizeRows(const DenseBase<T> &mat) {
    double row_sum;
    for (int i = 0; i < mat.rows(); ++i) {
        row_sum = mat.row(i).sum();
        if (row_sum == 0) { return false; }
        const_cast< DenseBase<T> & >(mat).row(i) /= row_sum;
    }
    return true;
}
SWIGCODE(%template(normalizeRows) normalizeRows<MatrixMd>;)

template<typename T>
bool normalizeCols(const DenseBase<T> &mat) {
    double col_sum;
    for (int j = 0; j < mat.cols(); ++j) {
        col_sum = mat.col(j).sum();
        if (col_sum == 0) { return false; }
        const_cast< DenseBase<T> & >(mat).col(j) /= col_sum;
    }
    return true;
}
SWIGCODE(%template(normalizeCols) normalizeCols<MatrixMd>;)

SWIGCODE(%apply const MatrixBase<MatrixXd>& OUTPUT { const MatrixBase<MatrixXd>& X };)

/**
 * Return the Kronecker-product \f$A\otimes B\f$ of the matrices `A` and `B`. */
template< typename D1, typename D2 >
MatrixXd kron(const MatrixBase<D1>& A, const MatrixBase<D2>& B) {
    MatrixXd result(A.rows() * B.rows(), A.cols() * B.cols());
    for (long i = 0; i < A.rows(); ++i) {
        for (long j = 0; j < A.cols(); ++j) {
            result.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i,j) * B;
        }
    }
    return result;
}
SWIGCODE(%template(kron) kron<MatrixMd, MatrixMd>;)

/**
 * Return
 * \ifnot PY
 * in the output-argument `X`,
 * \endif
 * the ordinary least-squares (OLS) solution to the problem `B` * `X` = `M`, using a `method` from {"Cholesky", "LDLT", "QR" (default), "SVD", "JacobiSVD"}.
 *
 * The "Cholesky" method solves the normal equations using a Cholesky decomposition. This is the fastest method, but loses most precision and requires the problem to be overdetermined and `B` to have full rank.
 *
 * The "LDLT" method is essentially the same as "Cholesky", but uses a more robust Cholesky decomposition with pivoting that also avoids taking a square root. This method is recommended over "Cholesky" by Eigen3.
 *
 * The "QR" method uses a QR decomposition. This is slower than "Cholesky", but gives more precision. The marix `B` should have full rank.
 *
 * The "SVD" uses an SVD decomposition. This is the slowest, but gives best precision. Also, the matrix `B` does not need to have full rank, and in the case of an underdetermined problem, the least-squares solution with the smallest norm is returned.
 *
 * The "JacobiSVD" method is similar to the "SVD" method, but uses a different (slower, but potentially more accurate) svd algorithm.
 */
template< typename D, typename D1, typename D2 >
C1(void) PY1(MatrixXd)
solveOLS(C2(const MatrixBase<D> &X,) const MatrixBase<D1>&B, const MatrixBase<D2>& M, const std::string& method = "QR") {
    if (method == "Cholesky") {
        const_cast<MatrixBase<D> &>(X) = (B.transpose() * B).llt().solve(B.transpose() * M);
    } else if (method == "LDLT") {
        const_cast<MatrixBase<D> &>(X) = (B.transpose() * B).ldlt().solve(B.transpose() * M);
    } else if (method == "QR") {
        const_cast<MatrixBase<D> &>(X) = B.colPivHouseholderQr().solve(M);
    } else if (method == "SVD") {
        const_cast<MatrixBase<D> &>(X) = B.bdcSvd(ComputeThinU | ComputeThinV).solve(M);
    } else if (method == "JacobiSVD") {
        const_cast<MatrixBase<D> &>(X) = B.jacobiSvd(ComputeThinU | ComputeThinV).solve(M);
    } else { const_cast<MatrixBase<D> &>(X).resize(0,0); }
}
SWIGCODE(%template(solveOLS) solveOLS<MatrixXd, MatrixXd, MatrixXd>;)

/**
 * Return
 * \ifnot PY
 * in the output-argument `X`,
 * \endif
 * the ordinary least-squares (OLS) solution to the problem `X` * `A` = `M`, using a `method` from {"Cholesky", "LDLT", "QR" (default), "SVD", "JacobiSVD"}.
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
tsolveOLS(C2(const MatrixBase<D> &X,) const MatrixBase<D1>& A, const MatrixBase<D2>& M, const std::string& method = "QR") {
    solveOLS(X, A.transpose(), M.transpose(), method);
    const_cast<MatrixBase<D> &>(X).transposeInPlace();
    return;
}
SWIGCODE(%template(tsolveOLS) tsolveOLS<MatrixXd, MatrixXd, MatrixXd>;)

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
    if (transpose) { tsolveOLS(result, M, MatrixXd::Identity(I_size, I_size), method); }
    else {            solveOLS(result, M, MatrixXd::Identity(I_size, I_size), method); }
    return result;
}
SWIGCODE(%template(pinv) pinv<MatrixMd>;)

/**
 * Return
 * \ifnot PY
 * in the output-argument `X`,
 * \endif
 * the (element-wise) D(`W`)-weighted least-squares (WLS) solution to the problem `B` * `X` = `M`, using a ` method` from {"Cholesky", "LDLT", "QR" (default), "SVD", "JacobiSVD"}, with weights given by `sqrt_W` such that `W` = `sqrt_W` .* `sqrt_W` where `.*` denotes the element-wise product.
 *
 * This computes `X` that minimizes || `W` .* (`B` * `X` - `M`) ||_F, where `.*` denotes the element-wise product by reducing the problem to a set of OLS problems that are then solved according to the given `method` as detailed below (see also `solveOLS()`). Note that the weights in `W` must be strictly greater than zero.
 *
 * The "Cholesky" method solves the normal equations using a Cholesky decomposition. This is the fastest method, but loses most precision and requires the problem to be overdetermined and `B` to have full rank.
 *
 * The "LDLT" method is essentially the same as "Cholesky", but uses a more robust Cholesky decomposition with pivoting that also avoids taking a square root. This method is recommended over "Cholesky" by Eigen3.
 *
 * The "QR" method uses a QR decomposition. This is slower than "Cholesky", but gives more precision. The marix `B` should have full rank.
 *
 * The "SVD" uses an SVD decomposition. This is the slowest, but gives best precision. Also, the matrix `B` does not need to have full rank, and in the case of an underdetermined problem, the least-squares solution with the smallest norm is returned.
 *
 * The "JacobiSVD" method is similar to the "SVD" method, but uses a different (slower, but potentially more accurate) svd algorithm.
 */
template< typename D, typename D1, typename D2, typename D3>
C1(void) PY1(MatrixXd)
solveW2LS(C2(const MatrixBase<D> &X,) const MatrixBase<D1> &B, const MatrixBase<D2> &M, const MatrixBase<D3> &sqrt_W, const std::string &method = "QR") {
    const_cast<MatrixBase<D>&>(X).derived().resize(B.cols(), M.cols());
    for (int j = 0; j < X.cols(); ++j) {
        solveOLS(const_cast<MatrixBase<D>&>(X).col(j), sqrt_W.col(j).asDiagonal() * B, sqrt_W.col(j).asDiagonal() * M.col(j), method);
    }
}
SWIGCODE(%template(solveW2LS) solveW2LS<MatrixXd, MatrixMd, MatrixMd, MatrixMd>;)

/**
 * Return
 * \ifnot PY
 * in the output-argument `X`,
 * \endif
 * the (element-wise) D(`W`)-weighted least-squares (WLS) solution to the problem `X` * `A` = `M`, using a ` method` from {"Cholesky", "LDLT", "QR" (default), "SVD", "JacobiSVD"}, with weights given by `sqrt_W` such that `W` = `sqrt_W` .* `sqrt_W` where `.*` denotes the element-wise product.
 *
 * This computes `X` that minimizes || `W` .* (`X` * `A` - `M`) ||_F, where `.*` denotes the element-wise product by reducing the problem to a set of OLS problems that are then solved according to the given `method` as detailed below (see also `solveOLS()`). Note that the weights in `W` must be strictly greater than zero.
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
tsolveW2LS(C2(const MatrixBase<D> &X,) const MatrixBase<D1> &A, const MatrixBase<D2> &M, const MatrixBase<D3> &sqrt_W, const std::string &method = "QR") {
    solveW2LS(X, A.transpose(), M.transpose(), sqrt_W.transpose(), method);
    const_cast<MatrixBase<D> &>(X).transposeInPlace();
    return;
}
SWIGCODE(%template(tsolveW2LS) tsolveW2LS<MatrixXd, MatrixMd, MatrixMd, MatrixMd>;)

/**
 * Return
 * \ifnot PY
 * in the output-argument `X`,
 * \endif
 * the (element-wise) D(`W`)-weighted least-squares (WLS) solution to the problem `B` * `X` = `M`, using a ` method` from {"Cholesky", "LDLT" (default), "QR", "SVD", "JacobiSVD"}.
 *
 * This computes `X` that minimizes || `W` .* (`B` * `X` - `M`) ||_F, where `.*` denotes the element-wise product by reducing the problem to a set of OLS problems that are then solved according to the given `method` as detailed below (see also `solveOLS()`). Note that the weights in `W` must be strictly greater than zero.
 *
 * The "Cholesky" method solves the normal equations using a Cholesky decomposition. This is the fastest method, but loses most precision and requires the problem to be overdetermined and `B` to have full rank.
 *
 * The "LDLT" method is essentially the same as "Cholesky", but uses a more robust Cholesky decomposition with pivoting that also avoids taking a square root. This method is recommended over "Cholesky" by Eigen3.
 *
 * The "QR" method uses a QR decomposition. This is slower than "Cholesky", but gives more precision. The marix `B` should have full rank.
 *
 * The "SVD" uses an SVD decomposition. This is the slowest, but gives best precision. Also, the matrix `B` does not need to have full rank, and in the case of an underdetermined problem, the least-squares solution with the smallest norm is returned.
 *
 * The "JacobiSVD" method is similar to the "SVD" method, but uses a different (slower, but potentially more accurate) svd algorithm.
 */
template< typename D, typename D1, typename D2, typename D3>
C1(void) PY1(MatrixXd)
solveWLS(C2(const MatrixBase<D> &X,) const MatrixBase<D1>& A, const MatrixBase<D2>& M, const MatrixBase<D3>& W, const std::string &method = "LDLT") {
    if (method == "Cholesky") {
        LLT<MatrixXd> llt;
        MatrixXd AT_WM = A.transpose() * (M.array() * W.array()).matrix();
        const_cast<MatrixBase<D>&>(X).derived().resize(A.cols(), M.cols());
        for (int j = 0; j < X.cols(); ++j) {
            llt.compute(A.transpose() * (W.col(j).asDiagonal() * A));
            const_cast<MatrixBase<D>&>(X).col(j).noalias() = llt.solve(AT_WM.col(j));
        }
    } else if (method == "LDLT") {
        LDLT<MatrixXd> ldlt;
        MatrixXd AT_WM = A.transpose() * (M.array() * W.array()).matrix();
        const_cast<MatrixBase<D>&>(X).derived().resize(A.cols(), M.cols());
        for (int j = 0; j < X.cols(); ++j) {
            ldlt.compute(A.transpose() * (W.col(j).asDiagonal() * A));
            const_cast<MatrixBase<D>&>(X).col(j).noalias() = ldlt.solve(AT_WM.col(j));
        }
    } else {
        solveW2LS(X, A, M, W.cwiseSqrt(), method);
    }
}
SWIGCODE(%template(solveWLS) solveWLS<MatrixXd, MatrixXd, MatrixXd, MatrixXd>;)

/**
 * Return
 * \ifnot PY
 * in the output-argument `X`,
 * \endif
 * the (element-wise) D(`W`)-weighted least-squares (WLS) solution to the problem `X` * `A` = `M`, using a ` method` from {"Cholesky", "LDLT" (default), "QR", "SVD", "JacobiSVD"}.
 *
 * This computes `X` that minimizes || `W` .* (`X` * `A` - `M`) ||_F, where `.*` denotes the element-wise product by reducing the problem to a set of OLS problems that are then solved according to the given `method` as detailed below (see also `solveOLS()`). Note that the weights in `W` must be strictly greater than zero.
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
tsolveWLS(C2(const MatrixBase<D> &X,) const MatrixBase<D1>& A, const MatrixBase<D2>& M, const MatrixBase<D3>& W, const std::string &method = "LDLT") {
    solveWLS(X, A.transpose(), M.transpose(), W.transpose(), method);
    const_cast<MatrixBase<D> &>(X).transposeInPlace();
    return;
}
SWIGCODE(%template(tsolveWLS) tsolveWLS<MatrixXd, MatrixXd, MatrixXd, MatrixXd>;)

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
            newWeights.middleCols(j * W.rows(), W.rows()).noalias() = B.transpose() * (W.col(j).asDiagonal() * B);
        } else {
            newWeights.col(j).noalias() = (B.transpose() * (W.col(j).asDiagonal() * B)).diagonal();
        }
    }
    return newWeights;
}

/**
 * Return
 * \ifnot PY
 * in the output-argument `X`,
 * \endif
 * the D(W1,..., Wm)-weighted least-squares (GLS) solution to the overdetermined problem `B` * `X` = `M`, using a `method` from {"Cholesky", "LDLT" (default)}, where the block-diagonal symmetric and positive definite weight matrix is given by `W` = [W1,..., Wm], where each `Wj` is the full weight matrix for the column j of `M`.
 *
 * This computes `X` that minimizes || S * (`B` * `X` - `M`) ||_F, where D(W1,..., Wm) = S^T * S is the Cholesky decomposition of the weight matrix.
 *
 * Note that the "LDLT" method is essentially the same as "Cholesky", but uses a more robust Cholesky decomposition with pivoting that also avoids taking a square root. This method is recommended over "Cholesky" by Eigen3.
 */
template< typename D, typename D1, typename D2, typename D3>
C1(void) PY1(MatrixXd)
solveGLS(C2(const MatrixBase<D> &X,) const MatrixBase<D1> &B, const MatrixBase<D2>& M, const MatrixBase<D3>& W, const std::string &method = "LDLT") {
    assert(W.rows() == M.rows() and W.cols() == M.cols() * M.rows()); // Block-diagonal weight matrix
    // solve AX≈M by solving SjAXj ≈ SjMj for each column j
    if (method == "Cholesky") {
        LLT<MatrixXd> llt;
        const_cast<MatrixBase<D>&>(X).derived().resize(B.cols(), M.cols());
        for (int j = 0; j < X.cols(); ++j) {
            llt.compute(B.transpose() * W.middleCols(j * W.rows(), W.rows()) * B);
            const_cast<MatrixBase<D>&>(X).col(j).noalias() = llt.solve(
                    B.transpose() * (W.middleCols(j * W.rows(), W.rows()) * M.col(j)));
        }
    } else if (method == "LDLT") {
        LDLT<MatrixXd> ldlt;
        const_cast<MatrixBase<D>&>(X).derived().resize(B.cols(), M.cols());
        for (int j = 0; j < X.cols(); ++j) {
            ldlt.compute(B.transpose() * W.middleCols(j * W.rows(), W.rows()) * B);
            const_cast<MatrixBase<D>&>(X).col(j).noalias() = ldlt.solve(
                    B.transpose() * (W.middleCols(j * W.rows(), W.rows()) * M.col(j)));
        }
    }
}
SWIGCODE(%template(solveGLS) solveGLS<MatrixXd, MatrixXd, MatrixXd, MatrixXd>;)

/**
 * Return
 * \ifnot PY
 * in the output-argument `X`,
 * \endif
 * the D(W1,..., Wm)-weighted least-squares (GLS) solution to the overdetermined problem `X` * `A` = `M`, using a `method` from {"Cholesky", "LDLT" (default)}, where the block-diagonal symmetric and positive definite weight matrix is given by `W` = [W1,..., Wm], where each `Wj` is the full weight matrix for the column j of `M`.
 *
 * This computes `X` that minimizes || S * (`X` * `A` - `M`) ||_F, where D(W1,..., Wm) = S^T * S is the Cholesky decomposition of the weight matrix.
 *
 * Note that the "LDLT" method is essentially the same as "Cholesky", but uses a more robust Cholesky decomposition with pivoting that also avoids taking a square root. This method is recommended over "Cholesky" by Eigen3.
 */
template< typename D, typename D1, typename D2, typename D3>
C1(void) PY1(MatrixXd)
tsolveGLS(C2(const MatrixBase<D> &X,) const MatrixBase<D1>& A, const MatrixBase<D2>& M, const MatrixBase<D3>& W, const std::string &method = "LDLT") {
    assert(W.rows() == M.rows() and W.cols() == M.cols() * M.rows()); // Block-diagonal weight matrix
    // solve XA≈M by solving [\sum_j AjAj^T ⊗ Wj] vec(X) = vec( [W1M1, W2M2, ..., MnMn] A^T )
    MatrixXd AtI_W_ATtI = MatrixXd::Zero(A.rows() * W.rows(), A.rows() * W.rows()); // we only compute lower triagonal part
    MatrixXd AiijAjjjWj = MatrixXd::Zero(W.rows(), W.rows());
    for (long jj = 0; jj < A.rows(); ++jj) {
        for (long ii = jj; ii < A.rows(); ++ii) {
            AiijAjjjWj.setZero(W.rows(), W.rows());
            for (long j = 0; j < M.cols(); ++j) {
                AiijAjjjWj += A(ii, j) * A(jj, j) * W.middleCols(j * W.rows(), W.rows());
            }
            AtI_W_ATtI.block(ii * W.rows(), jj * W.rows(), W.rows(), W.rows()) = AiijAjjjWj;
        }
    }
    MatrixXd WM = MatrixXd::Zero(W.rows(), M.cols());
    for (int j = 0; j < M.cols(); ++j) {
        WM.col(j).noalias() = W.middleCols(j*W.rows(), W.rows()) * M.col(j);
    }
    WM = WM * A.transpose();
    const Map<const VectorXd,1> vecWMAT(WM.data(),WM.size());
    if (method == "LDLT") { const_cast<MatrixBase<D>&>(X) = AtI_W_ATtI.ldlt().solve(vecWMAT); }
    else if (method == "Cholesky") { const_cast<MatrixBase<D>&>(X) = AtI_W_ATtI.llt().solve(vecWMAT); }
    else { return; }
    const_cast<MatrixBase<D>&>(X).derived().resize(M.rows(),A.rows());
}
SWIGCODE(%template(tsolveGLS) tsolveGLS<MatrixXd, MatrixXd, MatrixXd, MatrixXd>;)


/** return in the output arguments \c B and \c A (an approximation to) the best weighted rank-d approximation to \c M with element-wise weights \c W such that \f$\|BA - M\|_W\f$ is minimized.*/
template<typename D1, typename D2, typename D3, typename D4>
double improveWLRA(const MatrixBase<D1> &B, const MatrixBase<D2> &A, const MatrixBase<D3> &M,
                   const MatrixBase<D4> &W, double convergenceThreshold = 1e-5, int maxIterations = 100,
                   const std::string &method = "LDLT") {
    MatrixXd sqrtW = W.cwiseSqrt();
    int it = 0;
    double new_wrmse = std::numeric_limits<double>::infinity();
    double wrmse;
    do {
        wrmse = new_wrmse;
        if (method == "Cholesky" or method == "LDLT") {
            tsolveWLS(B, A, M, W, method);
            solveWLS(A, B, M, W, method);
        } else {
            tsolveW2LS(B, A, M, sqrtW, method);
            solveW2LS(A, B, M, sqrtW, method);
        }
        new_wrmse = (M - B*A).cwiseProduct(sqrtW).norm() / sqrt(M.size());
        it++;
    } while ( (wrmse - new_wrmse) > convergenceThreshold and it < maxIterations );
    return new_wrmse;
}
SWIGCODE(%template(improveWLRA) improveWLRA<MatrixMd, MatrixMd, MatrixMd, MatrixMd >;)

SWIGCODE(%clear const MatrixBase<MatrixXd>& X;)

} // namespace tom
