#include "Fun_assertion.h";

DyMat<double> DyMatMul(DyMat<double> A, DyMat<double> B) {
	if (A.yelm != B.xelm) {
		throw "the Dimension of the two matrix are not corresponding";
	}
	DyMat<double> D(A.xelm,B.yelm,(double)0);
	for (int i = 0; i < A.xelm; i++) {
		for (int j = 0; j < B.yelm; j++) {
			double T = 0;
			// caculate every element of the Matrix multiply
			for (int k = 0; k < B.xelm; k++) {
				T += A.get(i, k) * B.get(k, j);
			}
			D.set(i, j, T);
		}
	}
	return D;
};

/// calculate the inverse of the matrix by the Gaussian Eliminate Method
DyMat<double> DyMat_Inv(DyMat<double> A)
{
	if (A.xelm != A.yelm) {
		throw "The Matrix must be a square matrix";
	}
	int n = A.xelm;
	// Use the Gaussian elimination for the inverse of the matrix
	DyMat<double> W(n, 2 * n , (double)0 );
	DyMat<double> M(n , n, (double)0 );  // storge the result

	for (int i = 0; i < n ; i++) {
		for (int j = 0; j < n; j++) {
			// then we can extend the matrix
			W.set(i, j, A.get(i,j));
		}
		for (int j = n; j <2 * n; j++) {
			W.set(i, j, j - i == n ? 1 : 0);
		}
	}
	// note that the A.xelm and A.yelm has no difference 
	// firstly we have to judge if we can derive the inverse of the matrix
	for (int i = 0; i < n-1; i++) {  // the first line (original line)
		for (int j = n - 1; j > i; j--) {  // the second line (line need to change)
			if (abs(W.get(i, i)) < 0.00001) {
				throw "this matrix have no Inverse matrix!";
			}
			else {
				double coef = W.get(j, i) / W.get(i, i);
				// calculate the coefficient
				for (int k = 0; k < 2 * n; k++) {  // k is th column
					if (k <= i) {
						W.set(j, k, 0);
					}
					else {
						W.set(j, k, W.get(j, k) - W.get(i, k) * coef);
					}
					// transform every line of it as a unit matrix;
				}
			}
		}
	}
	// note that the condition that one matrix has its inverse matrix is that
	// it is a full-rank matrix ---> it can be derived by transforming it into a unit matrix
	if ( abs(W.get(n-1,n-1) < 0.00001) ) {
		cout << "Inverse_Matrix Error!";
		throw "the matrix have no inverse matrix";
	}
	else { // the matrix has inverse matrix, and calculate it;
		for (int i = 0; i < n-1; i++) {
			for (int j = i+1; j < n; j++) {
				double coef = W.get(i, j) / W.get(j, j);
				// i is the first line and the j is the second
				for (int k = j; k < 2 * n; k++) { // k is the column
					// I change k = j to k=n to save time
					if (k == j) {
						W.set(i, k, 0);  // in that case we will eliminate it directly
					}
					else {
						W.set(i, k, W.get(i, k) - W.get(j, k) * coef);
					}
				}
			}
		}
	}
	// W.Display_Content();
	// we can unify the previous matrix after that: 
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			M.set(i, j, W.get(i, j+n) / W.get(i, i));
		}
	}
	return M;
}

// DyMat<double> DyMatInv(DyMat<double> A,DyMat)
DyMat<double> DyMat_Add(DyMat<double> A, DyMat<double> B) {
	if (A.xelm != B.xelm || A.yelm != B.yelm) {
		throw "The dimension of two DyMatrix should be correspond!";
	}
	DyMat<double> C(A.xelm, A.yelm, (double)0);
	for (int i = 0; i < A.xelm; i++) {
		for (int j = 0; j < A.yelm; j++) {
			C.set(i, j, A.get(i, j) + B.get(i, j));
		}
	}
	return C;
}

DyMat<double> DyMat_Sub(DyMat<double> A, DyMat<double> B) {
	if (A.xelm != B.xelm || A.yelm != B.yelm) {
		throw "The dimension of two DyMatrix should be correspond!";
	}
	DyMat<double> C(A.xelm, A.yelm, (double)0);
	for (int i = 0; i < A.xelm; i++) {
		for (int j = 0; j < A.yelm; j++) {
			C.set(i, j, A.get(i, j) - B.get(i, j));
		}
	}
	return C;
}

DyMat<double> DyMat_dotMul(DyMat<double> A, DyMat<double> B) {
	if (A.xelm != B.xelm || A.yelm != B.yelm) {
		throw "The dimension of two DyMatrix should be correspond!";
	}
	DyMat<double> C(A.xelm, A.yelm, (double)0);
	for (int i = 0; i < A.xelm; i++) {
		for (int j = 0; j < A.yelm; j++) {
			C.set(i, j, A.get(i, j) * B.get(i, j));
		}
	}
	return C;
}

DyMat<double> DyMat_dotMul(DyMat<double> A, double b) {
	DyMat<double> C(A.xelm, A.yelm, (double)0);
	for (int i = 0; i < A.xelm; i++) {
		for (int j = 0; j < A.yelm; j++) {
			C.set(i, j, A.get(i, j) * b);
		}
	}
	return C;
}

DyMat<double> DyMat_dotDiv(DyMat<double> A, DyMat<double> B) {
	if (A.xelm != B.xelm || A.yelm != B.yelm) {
		throw "The dimension of two DyMatrix should be correspond!";
	}
	DyMat<double> C(A.xelm, A.yelm, (double)0);
	for (int i = 0; i < A.xelm; i++) {
		for (int j = 0; j < A.yelm; j++) {
			C.set(i, j, A.get(i, j) / B.get(i, j));
		}
	}
	return C;
}

DyMat<double> DyMat_dotDiv(DyMat<double> A, double b) {
	DyMat<double> C(A.xelm, A.yelm, (double)0);
	for (int i = 0; i < A.xelm; i++) {
		for (int j = 0; j < A.yelm; j++) {
			C.set(i, j, A.get(i, j) / b);
		}
	}
	return C;
}

DyMat<double> DyMat_dotPow(DyMat<double> A, double b) {
	DyMat<double> C(A.xelm, A.yelm, (double)0);
	for (int i = 0; i < A.xelm; i++) {
		for (int j = 0; j < A.yelm; j++) {
			C.set(i, j, pow(A.get(i, j),b));
		}
	}
	return C;
}

DyMat<double> DyMat_transpose(DyMat<double> A) {
	DyMat<double> C(A.yelm, A.xelm, (double)0);
	for (int i = 0; i < A.xelm; i++) {
		for (int j = 0; j < A.yelm; j++) {
			C.set(j, i, A.get(i, j));
		}
	}
	return C;
}

double DyMat_sum(DyMat<double> A) {
	double _sum=0;
	for (int i = 0; i < A.xelm; i++) {
		for (int j = 0; j < A.yelm; j++) {
			_sum += A.get(i, j);
		}
	}
	return _sum;
}