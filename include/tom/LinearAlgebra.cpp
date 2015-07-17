#include "tom.h"

namespace tom {

#ifndef SWIG
typedef Map<Matrix< double, Dynamic, Dynamic>, 0, Stride<Dynamic, Dynamic> >  MatrixMd;
#endif

template< typename D1, typename D2 >
void pinvFast(const Eigen::MatrixBase<D1>& result, const Eigen::MatrixBase<D2>& M, const std::string& method) {
	if (M.rows() >= M.cols()) {
		if (method == "LDLT") {
			const_cast<MatrixBase<D1>&>(result) = (M.transpose() * M).ldlt().solve(M.transpose());
		} else { // method LLT as default
			const_cast<MatrixBase<D1>&>(result) = (M.transpose() * M).llt().solve(M.transpose());
		}
	} else { // M.rows() < M.cols()
		if (method == "LDLT") {
			const_cast<MatrixBase<D1>&>(result) = (M * M.transpose()).ldlt().solve(M).transpose();
		} else { // method LLT as default
			const_cast<MatrixBase<D1>&>(result) = (M * M.transpose()).llt().solve(M).transpose();
		}
	}
}
// template void pinvFast(const MatrixBase<MatrixXd>&, const MatrixBase<MatrixMd>&, const std::string& );
// template void pinvFast(const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, const std::string& );


template< typename D1, typename D2 >
void pinv(const MatrixBase<D1>& result, const MatrixBase<D2>& M, double tolerance, const std::string& method) {
	if (method == "SVD") {
		JacobiSVD<MatrixXd> svd(M, ComputeThinU | ComputeThinV);
		VectorXd s = svd.singularValues();
		if (tolerance == -1)
			tolerance = std::max(M.rows(), M.cols()) * s(0) * NumTraits<double>::epsilon();
		long rank = 0;
		for ( long i = 0; i < s.size(); ++i)
			if (s.coeff(i) > tolerance) { rank++; s.coeffRef(i) = 1.0 / s.coeff(i); }
		const_cast<MatrixBase<D1>&>(result).noalias() = svd.matrixV().leftCols(rank) * s.head(rank).asDiagonal() * svd.matrixU().leftCols(rank).transpose();
	}
	else { // method BDCSVD as default
		BDCSVD<MatrixXd> svd(M, ComputeThinU | ComputeThinV);
		VectorXd s = svd.singularValues();
		if (tolerance == -1)
			tolerance = std::max(M.rows(), M.cols()) * s(0) * NumTraits<double>::epsilon();
		long rank = 0;
		for ( long i = 0; i < s.size(); ++i)
			if (s.coeff(i) > tolerance) { rank++; s.coeffRef(i) = 1.0 / s.coeff(i); }
		const_cast<MatrixBase<D1>&>(result).noalias() = svd.matrixV().leftCols(rank) * s.head(rank).asDiagonal() * svd.matrixU().leftCols(rank).transpose();
	}
}
//template void pinv(const MatrixBase<MatrixXd>&, const MatrixBase<MatrixMd>&, double, const std::string& );
//template void pinv(const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, double, const std::string& );


template< typename D >
MatrixXd pinv(const MatrixBase<D>& M, double tolerance) {
	MatrixXd result;
	pinv(result, M, tolerance);
	return result;
}
//template MatrixXd pinv(const MatrixBase<MatrixXd>&, double);


void eigensolve(VectorXcd& eigenvals, MatrixXcd& eigenvecs, const MatrixXd& M) {
	EigenSolver<MatrixXd> eigensolver(M);
	eigenvals = eigensolver.eigenvalues();
	eigenvecs = eigensolver.eigenvectors();
}


template< typename D1, typename D2, typename D3 >
void solveFastOLS(const MatrixBase<D1>& X, const MatrixBase<D2>& A, const MatrixBase<D3>& M, bool transposed, const std::string& method) {
	if (method == "LDLT") {
		if (transposed)
			const_cast<MatrixBase<D1>&>(X) = (A.transpose() * A).ldlt().solve(A.transpose() * M);
		else
			const_cast<MatrixBase<D1>&>(X) = (A * A.transpose()).ldlt().solve(A * M.transpose()).transpose();
	} else { // method LLT as default
		if (transposed)
			const_cast<MatrixBase<D1>&>(X) = (A.transpose() * A).llt().solve(A.transpose() * M);
		else
			const_cast<MatrixBase<D1>&>(X) = (A * A.transpose()).llt().solve(A * M.transpose()).transpose();
	}
}
//template void solveFastOLS(const MatrixBase<MatrixXd>&, const MatrixBase<MatrixMd>&, const MatrixBase<MatrixMd>&, bool, const std::string&);
//template void solveFastOLS(const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, bool, const std::string&);


template< typename D1, typename D2, typename D3 >
void solveOLS(const MatrixBase<D1>& X, const MatrixBase<D2>& A, const MatrixBase<D3>& M, bool transposed, const std::string& method) {
	if (method == "LDLT") {
		if (transposed)
			const_cast<MatrixBase<D1>&>(X) = (A.transpose() * A).ldlt().solve(A.transpose() * M);
		else
			const_cast<MatrixBase<D1>&>(X) = (A * A.transpose()).ldlt().solve(A * M.transpose()).transpose();
	} else if (method == "LLT") {
		if (transposed)
			const_cast<MatrixBase<D1>&>(X) = (A.transpose() * A).llt().solve(A.transpose() * M);
		else
			const_cast<MatrixBase<D1>&>(X) = (A * A.transpose()).llt().solve(A * M.transpose()).transpose();
	}	else if (method == "SVD") {
		if (transposed)
			const_cast<MatrixBase<D1>&>(X) = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(M);
		else
			const_cast<MatrixBase<D1>&>(X) = A.transpose().jacobiSvd(ComputeThinU | ComputeThinV).solve(M.transpose()).transpose();

	} else if (method == "QR") {
		if (transposed)
			const_cast<MatrixBase<D1>&>(X) = A.colPivHouseholderQr().solve(M);
		else
			const_cast<MatrixBase<D1>&>(X) = A.transpose().colPivHouseholderQr().solve(M.transpose()).transpose();
	} else { // method BDCSVD as default
		if (transposed)
			const_cast<MatrixBase<D1>&>(X) = A.bdcSvd(ComputeThinU | ComputeThinV).solve(M);
		else
			const_cast<MatrixBase<D1>&>(X) = A.transpose().bdcSvd(ComputeThinU | ComputeThinV).solve(M.transpose()).transpose();
	}
}
//template void solveOLS(const MatrixBase<MatrixXd>&, const MatrixBase<MatrixMd>&, const MatrixBase<MatrixMd>&, bool, const std::string&);
//template void solveOLS(const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, bool, const std::string&);


template< typename D1, typename D2 , typename D3, typename D4>
void solveFastWLS(const MatrixBase<D1>& X, const MatrixBase<D2>& A, const MatrixBase<D3>& M, const MatrixBase<D4>& W, bool transposed, const std::string& method) {
	if (transposed) { const_cast<MatrixBase<D1>&>(X).derived().resize(A.cols(), M.cols()); }
	else { const_cast<MatrixBase<D1>&>(X).derived().resize(M.rows(), A.rows()); }
	if (method == "LDLT") {
		if (transposed) {
			for (int j = 0; j < X.cols(); ++j) {
				const_cast<MatrixBase<D1>&>(X).col(j) = (A.transpose() * W.col(j).asDiagonal() * A).ldlt().solve(A.transpose() * W.col(j).asDiagonal() * M.col(j));
			}
		} else {
			for (int i = 0; i < X.rows(); ++i) {
				const_cast<MatrixBase<D1>&>(X).row(i) = (A * W.row(i).asDiagonal() * A.transpose()).ldlt().solve(A * W.row(i).asDiagonal() * M.row(i).transpose()).transpose();
			}
		}
	} else { // method LLT as defalut
		if (transposed) {
			for (int j = 0; j < X.cols(); ++j) {
				const_cast<MatrixBase<D1>&>(X).col(j) = (A.transpose() * W.col(j).asDiagonal() * A).llt().solve(A.transpose() * W.col(j).asDiagonal() * M.col(j));
			}
		} else {
			for (int i = 0; i < X.rows(); ++i) {
				const_cast<MatrixBase<D1>&>(X).row(i) = (A * W.row(i).asDiagonal() * A.transpose()).llt().solve(A * W.row(i).asDiagonal() * M.row(i).transpose()).transpose();
			}
		}
	}
}
//template void solveFastWLS(const MatrixBase<MatrixXd>&, const MatrixBase<MatrixMd>&, const MatrixBase<MatrixMd>&, const MatrixBase<MatrixMd>&, bool, const std::string&);
//template void solveFastWLS(const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, bool, const std::string&);


template< typename D1, typename D2 , typename D3, typename D4>
void solveWLS(const MatrixBase<D1>& X, const MatrixBase<D2>& A, const MatrixBase<D3>& M, const MatrixBase<D4>& sqrtW, bool transposed, const std::string& method) {
	if (transposed) { const_cast<MatrixBase<D1>&>(X).derived().resize(A.cols(), M.cols()); }
	else { const_cast<MatrixBase<D1>&>(X).derived().resize(M.rows(), A.rows()); }
	if (transposed) {
		for (int j = 0; j < X.cols(); ++j) {
			solveOLS(const_cast<MatrixBase<D1>&>(X).col(j), sqrtW.col(j).asDiagonal() * A, sqrtW.col(j).asDiagonal() * M.col(j), transposed, method);
		}
	} else {
		for (int i = 0; i < X.rows(); ++i) {
			solveOLS(const_cast<MatrixBase<D1>&>(X).row(i), A * sqrtW.row(i).asDiagonal(), M.row(i) * sqrtW.row(i).asDiagonal(), transposed, method);
		}
	}
}
//template void solveWLS(const MatrixBase<MatrixXd>&, const MatrixBase<MatrixMd>&, const MatrixBase<MatrixMd>&, const MatrixBase<MatrixMd>&, bool, const std::string&);
//template void solveWLS(const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, bool, const std::string&);


template< typename D1, typename D2 >
void solveFastGLS(Matrix<double, Dynamic, Dynamic, ColMajor>& X,
									const MatrixBase<D1>& A, const Matrix<double, Dynamic, Dynamic, ColMajor>& M, const MatrixBase<D2>& W, bool transposed, const std::string& method) {
	if (W.rows() == M.rows()) { // block-diagonal weight matrix
		int m = M.rows();
		if (transposed) { // solve AX≈M by solving SjAXj ≈ SjMj for each column j
			X.resize(A.cols(), M.cols());
			if (method == "LDLT") {
				for (int j = 0; j < X.cols(); ++j) {
					X.col(j) = (A.transpose() * W.middleCols(j*m, m) * A).ldlt().solve(A.transpose() * W.middleCols(j*m, m) * M.col(j));
				}
			} else { // defalut method LLT
				for (int j = 0; j < X.cols(); ++j) {
					X.col(j) = (A.transpose() * W.middleCols(j*m, m) * A).llt().solve(A.transpose() * W.middleCols(j*m, m) * M.col(j));
				}
			}
		} else { // solve XA≈M by solving [\sum_j AjAj^T ⊗ Wj] vec(X) = vec( [W1M1, W2M2, ..., MnMn] A^T )
			MatrixXd AtI_W_ATtI = MatrixXd::Zero(A.rows()*m, A.rows()*m);
			MatrixXd WM = MatrixXd::Zero(M.rows(), M.cols());
			for (int j = 0; j < M.cols(); ++j) {
				AtI_W_ATtI += kron(A.col(j) * A.col(j).transpose(), W.middleCols(j*m, m));
				WM.col(j) = W.middleCols(j*m, m) * M.col(j);
			}
			WM *= A.transpose();
			const Map<const VectorXd> vecWMAT(WM.data(),WM.size());
			if (method == "LDLT") { X = AtI_W_ATtI.ldlt().solve(vecWMAT); }
			else { /* default LLT */ X = AtI_W_ATtI.llt().solve(vecWMAT); }
			X.resize(m,A.rows());
		}
	} else { // full weight matrix
		const Map<const VectorXd> vecM(M.data(),M.size());
		if (transposed) { // solve AX≈M by solving S(I ⊗ A)vec(X) ≈ Svec(M)
			MatrixXd ItA = kron(MatrixXd::Identity(M.cols(), M.cols()), A);
			if (method == "LDLT") {
				const_cast<Matrix<double, Dynamic, Dynamic, ColMajor>&>(X) = (ItA.transpose() * W * ItA).ldlt().solve(ItA.transpose() * W * vecM);
			} else if (method == "LLT") {
				const_cast<Matrix<double, Dynamic, Dynamic, ColMajor>&>(X) = (ItA.transpose() * W * ItA).llt().solve(ItA.transpose() * W * vecM);
			}
			X.resize(A.cols(), M.cols());
		} else { // solve XA≈M by solving S(A^T ⊗ I)vec(X) ≈ Svec(M)
			MatrixXd ATtI = kron(A.transpose(), MatrixXd::Identity(M.rows(), M.rows()));
			if (method == "LDLT") {
				const_cast<Matrix<double, Dynamic, Dynamic, ColMajor>&>(X) = (ATtI.transpose() * W * ATtI).ldlt().solve(ATtI.transpose() * W * vecM);
			} else if (method == "LLT") {
				const_cast<Matrix<double, Dynamic, Dynamic, ColMajor>&>(X) = (ATtI.transpose() * W * ATtI).llt().solve(ATtI.transpose() * W * vecM);
			}
			X.resize(M.rows(), A.rows());
		}
	}
}
//template void solveFastGLS(Matrix<double, Dynamic, Dynamic, ColMajor>&, const MatrixBase<MatrixMd>&, const Matrix<double, Dynamic, Dynamic, ColMajor>&, const MatrixBase<MatrixMd>&, bool, const std::string&);
//template void solveFastGLS(Matrix<double, Dynamic, Dynamic, ColMajor>&, const MatrixBase<MatrixXd>&, const Matrix<double, Dynamic, Dynamic, ColMajor>&, const MatrixBase<MatrixXd>&, bool, const std::string&);



template< typename D1, typename D2 >
void solveGLS(Matrix<double, Dynamic, Dynamic, ColMajor>& X,
							const MatrixBase<D1>& A, const Matrix<double, Dynamic, Dynamic, ColMajor>& M, const MatrixBase<D2>& S, bool transposed, const std::string& method) {
	const Map<const VectorXd> vecM(M.data(),M.size());
	if (transposed) { // solve AX≈M by solving S(I ⊗ A)vec(X) ≈ Svec(M)
		MatrixXd ItA; kron(ItA, MatrixXd::Identity(M.cols(), M.cols()), A);
		solveOLS(X, S * ItA, S * vecM, true, method);
		X.resize(A.cols(), M.cols());
	} else { // solve XA≈M by solving S(A^T ⊗ I)vec(X) ≈ Svec(M)
		MatrixXd ATtI; kron(ATtI, A.transpose(), MatrixXd::Identity(M.rows(), M.rows()));
		solveOLS(X, S * ATtI, S * vecM, true, method);
		X.resize(M.rows(), A.rows());
	}
}
//template void solveGLS(Matrix<double, Dynamic, Dynamic, ColMajor>&, const MatrixBase<MatrixMd>&, const Matrix<double, Dynamic, Dynamic, ColMajor>&, const MatrixBase<MatrixMd>&, bool, const std::string&);
//template void solveGLS(Matrix<double, Dynamic, Dynamic, ColMajor>&, const MatrixBase<MatrixXd>&, const Matrix<double, Dynamic, Dynamic, ColMajor>&, const MatrixBase<MatrixXd>&, bool, const std::string&);


template< typename D1, typename D2 , typename D3, typename D4>
double improveWLRA(const MatrixBase<D1>& B, const MatrixBase<D2>& A, const MatrixBase<D3>& M, const MatrixBase<D4>& W, double convergenceThreshold, int maxIterations, const std::string& method) {
	MatrixXd sqrtW = W.cwiseSqrt();
	int it = 0;
	double new_wrmse = std::numeric_limits<double>::infinity();
	double wrmse;
	do {
		wrmse = new_wrmse;
		if (method == "LLT" or method == "LDLT") {
			solveFastWLS(B, A, M, W, false, method);
			solveFastWLS(A, B, M, W, true, method);
		} else {
			solveWLS(B, A, M, sqrtW, false, method);
			solveWLS(A, B, M, sqrtW, true, method);
		}
		new_wrmse = (M - B*A).cwiseProduct(sqrtW).norm() / sqrt(M.size());
		it++;
	} while ( (wrmse - new_wrmse) > convergenceThreshold and it < maxIterations );
	return new_wrmse;
}
//template double improveWLRA(const MatrixBase<MatrixMd>&, const MatrixBase<MatrixMd>&, const MatrixBase<MatrixMd>&, const MatrixBase<MatrixMd>&, double, int, const std::string&);
//template double improveWLRA(const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, const MatrixBase<MatrixXd>&, double, int, const std::string&);

// template< typename D1, typename D2, typename D3, typename D4 >
// MatrixXd EW_TLS(const MatrixBase<D1>& A, const MatrixBase<D2>& B, const MatrixBase<D3>& WA, const MatrixBase<D4>& WB) {
// 	MatrixXd AB, W;
// 	AB.resize(A.rows() + B.rows(), A.cols()); W.resize(A.rows() + B.rows(), A.cols());
// 	AB.topRows(A.rows()) = A; AB.bottomRows(B.rows()) = B;
// 	W.topRows(A.rows()) = WA; W.bottomRows(B.rows()) = WB;
// 	MatrixXd Bret, Aret;
// 	weightedLowRankApproximation(Bret, Aret, AB, W, A.rows());
// 	return Bret.bottomRows(B.rows()) * pinv(Bret.topRows(A.rows()));
// }


} // namespace tom
