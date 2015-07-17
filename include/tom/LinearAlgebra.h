#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include "tom.h"

namespace tom {

/**
 * return \f$a^b\f$ for non-negative integer bases and exponents.\ Note that \f$0^0 := 1\f$.
 */
inline
int ipow(int a, int b) {
	int res = 1;
	while (b) {
		if (b & 1) res *= a;
		b >>= 1;
		a *= a;
	}
	return res;
}

template <typename T>
double normalize(const Eigen::DenseBase<T>& mat) {
	double mat_sum = mat.sum();
	if (mat_sum == 0) { return 0; }
	const_cast< Eigen::DenseBase<T>& >(mat) /= mat_sum;
	return mat_sum;
}
TEMPLATE(normalize, normalize, MatrixMd)

template <typename T>
bool normalizeRows(const Eigen::DenseBase<T>& mat) {
	double row_sum;
	for (int i = 0; i < mat.rows(); ++i) {
		row_sum = mat.row(i).sum();
		if (row_sum == 0) { return false; }
		const_cast< Eigen::DenseBase<T>& >(mat).row(i) /= row_sum;
	}
	return true;
}
TEMPLATE(normalizeRows, normalizeRows, MatrixMd)

template <typename T>
bool normalizeCols(const Eigen::DenseBase<T>& mat) {
	double col_sum;
	for (int j = 0; j < mat.cols(); ++j) {
		col_sum = mat.col(j).sum();
		if (col_sum == 0) { return false; }
		const_cast< Eigen::DenseBase<T>& >(mat).col(j) /= col_sum;
	}
	return true;
}
TEMPLATE(normalizeCols, normalizeCols, MatrixMd)

/**
 * return in the output argument \c AB the Kronecker-product \f$A\otimes B\f$ of the matrices \c A and \c B.
 */
OUTMAP(AB)
template <typename D1, typename D2, typename D3>
void kron(const Eigen::MatrixBase<D1>& AB, const Eigen::MatrixBase<D2>& A, const Eigen::MatrixBase<D3>& B) {
	const_cast<Eigen::MatrixBase<D1>&>(AB).derived().resize(A.rows() * B.rows(), A.cols() * B.cols());
	for (int i = 0; i < A.rows(); ++i) {
		for (int j = 0; j < A.cols(); ++j) {
			const_cast<Eigen::MatrixBase<D1>&>(AB).block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i,j) * B;
		}
	}
}
TEMPLATE(kron, kron, MatrixXd, MatrixMd, MatrixMd)
CLEAROUTMAP(result)

template< typename D1, typename D2 >
Eigen::MatrixXd kron(const Eigen::MatrixBase<D1>& A, const Eigen::MatrixBase<D2>& B) {
	Eigen::MatrixXd result;
	kron(result, A, B);
	return result;
}

// template <typename Derived1, typename Derived2>
// void
// k_kron_times_vec_inplace(const Eigen::MatrixBase<Derived1>& A, const Eigen::MatrixBase<Derived2>& v, int k) {
// 	// The following is a standard Eigen hack to make v editable:
// 	Eigen::MatrixBase<Derived2>& v_ =	const_cast<Eigen::MatrixBase<Derived2>& >(v);
// 	const int n = A.cols();
// 	for (int it = 0; it < k; ++it) {
// 		int v_rows = ipow(n, it);
// 		int v_cols = ipow(n, k-it);
// 		int v_A_segs = ipow(n, k-it-1);
// 		Eigen::Map<Eigen::MatrixXd> v_mat(v_.derived().data(), v_rows, v_cols);
// 		for (int i = 0; i < v_A_segs; ++i) {
// 			v_mat.block(0, i*n, v_rows, n) = v_mat.block(0, i*n, v_rows, n) * (A.transpose());
// 		}
// 	}
// }
// #ifdef SWIG
// %template(k_kron_times_vec_inplace) k_kron_times_vec_inplace<MatrixMd, MatrixMd>;
// #endif


/**
 * for the given matrix \c M of size m x n, return in the output argument \c result the matrix \f$(M^\top M)^{-1}M^\top\f$ if m >= n, or \f$M^\top(MM^\top)^{-1}\f$ if m < n, which is computed as the OLS solution X to \f$MX\approx I\f$, or \f$XM\approx I\f$, respectively, using a \c method from {"LLT" (default), "LDLT"}\. \b Note that for this the matrix \c M must not be rank deficient!
 */
OUTMAP(result)
template< typename D1, typename D2 >
void pinvFast(const Eigen::MatrixBase<D1>& result, const Eigen::MatrixBase<D2>& M, const std::string& method = "LLT");
TEMPLATE(pinvFast, pinvFast, MatrixXd, MatrixMd)
CLEAROUTMAP(result)


/**
 * return in the output argument \c result the pseudoinverse \f$M^\dagger = \f$ of the given matrix \c M, computed using a \c method from {"BDCSVD" (default), "SVD"}\. All singular values below \c tolerance are treated as zero\. If \c tolerance = -1, the default tolerance is used.
 */
OUTMAP(result)
template< typename D1, typename D2 >
void pinv(const Eigen::MatrixBase<D1>& result, const Eigen::MatrixBase<D2>& M, double tolerance = -1, const std::string& method = "BDCSVD");
TEMPLATE(pinv, pinv, MatrixXd, MatrixMd)
CLEAROUTMAP(result)


template< typename D >
Eigen::MatrixXd pinv(const Eigen::MatrixBase<D>& M, double tolerance = -1);


/** return the complex eigenvalues and eigenvectors of the given matrix \c M in the output arguments \c eigvals and \c eigvecs. */
void eigensolve(Eigen::VectorXcd& eigenvals, Eigen::MatrixXcd& eigenvecs, const Eigen::MatrixXd& M);


/** return in the output argument \c X the OLS solution the overdetermined problem \f$XA\approx M\f$, or if \c transposed is set to \c true, to the overdetermined problem \f$AX\approx M\f$, using a \c method from {"LLT" (default), "LDLT"}. **/
OUTMAP(X)
template< typename D1, typename D2, typename D3 >
void solveFastOLS(const Eigen::MatrixBase<D1>& X, const Eigen::MatrixBase<D2>& A, const Eigen::MatrixBase<D3>& M, bool transposed = false, const std::string& method = "LLT");
TEMPLATE(solveFastOLS, solveFastOLS, MatrixXd, MatrixMd, MatrixMd)
CLEAROUTMAP(X)


/** return in the output argument \c X the OLS solution to the problem \f$XA\approx M\f$, or if \c transposed is set to \c true, to the problem \f$AX\approx M\f$, using a \c method from {"BDCSVD" (default), "SVD"}. **/
OUTMAP(X)
template< typename D1, typename D2, typename D3 >
void solveOLS(const Eigen::MatrixBase<D1>& X, const Eigen::MatrixBase<D2>& A, const Eigen::MatrixBase<D3>& M, bool transposed = false, const std::string& method = "BDCSVD");
TEMPLATE(solveOLS, solveOLS, MatrixXd, MatrixMd, MatrixMd)
CLEAROUTMAP(X)

/** return in the output argument \c X the WLS solution the overdetermined problem \f$XA\approx M\f$, or if \c transposed is set to \c true, to the overdetermined problem \f$AX\approx M\f$, with element-wise weights \c W, using a \c method from {"LLT" (default), "LDLT"}. **/
OUTMAP(X)
template< typename D1, typename D2 , typename D3, typename D4>
void solveFastWLS(const Eigen::MatrixBase<D1>& X, const Eigen::MatrixBase<D2>& A, const Eigen::MatrixBase<D3>& M, const Eigen::MatrixBase<D4>& W, bool transposed = false, const std::string& method = "LLT");
TEMPLATE(solveFastWLS, solveFastWLS, MatrixXd, MatrixMd, MatrixMd, MatrixMd)
CLEAROUTMAP(X)

/** return in the output argument \c X the WLS solution the problem \f$XA\approx M\f$, or if \c transposed is set to \c true, to the problem \f$AX\approx M\f$, with element-wise weights W given by \c sqrtW such that (\c sqrtW)^2 = W, using a \c method from {"BDCSVD" (default), "SVD"}. **/
OUTMAP(X)
template< typename D1, typename D2 , typename D3, typename D4>
void solveWLS(const Eigen::MatrixBase<D1>& X, const Eigen::MatrixBase<D2>& A, const Eigen::MatrixBase<D3>& M, const Eigen::MatrixBase<D4>& sqrtW, bool transposed = false, const std::string& method = "BDCSVD");
TEMPLATE(solveWLS, solveWLS, MatrixXd, MatrixMd, MatrixMd, MatrixMd)
CLEAROUTMAP(X)

/** return in the output argument \c X the GLS solution the overdetermined problem \f$XA\approx M\f$, or if \c transposed is set to \c true, to the overdetermined problem \f$AX\approx M\f$, with weights \c W, using a \c method from {"LLT" (default), "LDLT"}\. Note that \c W will in general be a matrix of size mn x mn for a matrix \c M of size m x n\. However, if a matrix \c W of size m x mn is passed, it is assumed that \c W = [W_1, W_2, ..., W_n] and the block-diagonal weight matrix D(W_1, W_2, ..., W_n) will be used. **/
OUTMATRIXXD(X)
template< typename D1, typename D2 >
void solveFastGLS(Eigen::MatrixXd& X, const Eigen::MatrixBase<D1>& A, const Eigen::MatrixXd& M, const Eigen::MatrixBase<D2>& W, bool transposed = false, const std::string& method = "LLT");
TEMPLATE(solveFastGLS, solveFastGLS, MatrixMd, MatrixMd)
CLEAROUTMATRIXXD(X)

/** return in the output argument \c X the GLS solution the overdetermined problem \f$XA\approx M\f$, or if \c transposed is set to \c true, to the overdetermined problem \f$AX\approx M\f$, with weights W given by \c S such that \f$S^\top S = W\f$, using a \c method from {"BDCSVD" (default), "SVD"}. **/
OUTMATRIXXD(X)
template< typename D1, typename D2 >
void solveGLS(Eigen::MatrixXd& X, const Eigen::MatrixBase<D1>& A, const Eigen::MatrixXd& M, const Eigen::MatrixBase<D2>& S, bool transposed = false, const std::string& method = "BDCSVD");
TEMPLATE(solveGLS, solveGLS, MatrixMd, MatrixMd)
CLEAROUTMATRIXXD(X)

/** return in the output arguments \c B and \c A (an approximation to) the best weighted rank-d approximation to \c M with element-wise weights \c W such that \f$\|BA - M\|_W\f$ is minimized.*/
template< typename D1, typename D2 , typename D3, typename D4>
double improveWLRA(const Eigen::MatrixBase<D1>& B, const Eigen::MatrixBase<D2>& A, const Eigen::MatrixBase<D3>& M, const Eigen::MatrixBase<D4>& W, double convergenceThreshold = 1e-5, int maxIterations = 100, const std::string& method = "BDCSVD");
TEMPLATE(improveWLRA, improveWLRA, MatrixMd, MatrixMd, MatrixMd, MatrixMd)


// /** return the the elementwise-weighted total least squares (EW-TLS) solution \c X for \f$XA \approx B\f$ with weights \c WA and \c WB.\ Note that it is assumed that \c A and \c B are of size n-times-m with n << m and that \c A has full row rank. */
// template< typename D1, typename D2, typename D3, typename D4 >
// Eigen::MatrixXd EW_TLS(const Eigen::MatrixBase<D1>& A, const Eigen::MatrixBase<D2>& B, const Eigen::MatrixBase<D3>& WA, const Eigen::MatrixBase<D4>& WB);

} // namespace tom

#endif // LINEAR_ALGEBRA_H
