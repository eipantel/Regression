#include "mmatrix.h"


using namespace std;

Matrix::Matrix(const Matrix& M) :
    elements(new double[M.n * M.m]),
    n(M.n),
    m(M.m)
{
    for (int i = 0; i < n*m; i++)
        elements[i] = M.elements[i];
}

Matrix::Matrix(int rows, int columns) :
    elements(new double[rows * columns]),
    n(rows),
    m(columns)
{
	for (int i = 0; i < n*m; i++)
		elements[i] = 0;
}

Matrix::Matrix(const int rows, const int columns, double* el) :
    elements(new double[rows * columns]),
    n(rows),
    m(columns)
{
    for (int i = 0; i < n*m; i++)
        elements[i] = el[i];
}

double& Matrix::at(int i, int j) {
    return elements[i*m + j];
}

double Matrix::at(int i, int j) const {				//??????????? ?????
	return elements[i*m + j];
}

void Matrix::operator=(const Matrix& mat) {    // Copy-operator
	for (int i = 0; i < m*n; i++) {
		elements[i] = mat.elements[i];
	}
}

Matrix Matrix::transpose() const {
	Matrix a(m, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			a.at(j, i) = at(i, j);
		}
	}
	return a;
}

Matrix Matrix::operator*(double k) const {
	Matrix ans(n, m);
	for (int i = 0; i < n*m; ++i) {
		ans.elements[i] = elements[i]*k;
	}
	return ans;
}

Matrix Matrix::operator/(double k) const {
    return (*this)*(1/k);
}

Matrix Matrix::operator+(const Matrix& A) const{
	if (A.n != n || A.m != m) throw "Matrix::operator+ - wrong dimention";
	Matrix ans(n, m);
	for (int i = 0; i < n*m; ++i) {
		ans.elements[i] = elements[i]+A.elements[i];
	}
	return ans;
}

Matrix Matrix::operator-(const Matrix& A) const {
	if (A.n != n || A.m != m) throw "Matrix::operator+ - wrong dimention";
	return *this + (A*(double)(-1));
}

Matrix Matrix::operator*(const Matrix& A) const{
	if (m != A.n) throw "Matrix::operator* - wrong dimention";
	Matrix ans(n, A.m);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < A.m; ++j) {
			for (int k = 0; k < m; ++k) ans.at(i, j) += (double)(at(i, k) * A.at(k, j));
		}
	}
	return ans;
}

void Matrix::ch_rows(int row1, int row2) {
    for (int i = 0; i < m; ++i) {
        double tmp = at(row2, i);
        at(row2, i) = at(row1, i);
        at(row1, i) = tmp;
    }
}

double Matrix::det() const {
    if (n != m) throw "Matrix::det - wrong dimension";
    Matrix mat(n, m, elements);
    if (n != m) return 0;
    else {
        double det = 1;
        for (int i = 0; i < n; ++i) {
            if (mat.at(i, i) == 0) {
                bool flag = true;
                for (int j = i+1; j < n; ++j) {
                    if (mat.at(j, i) != 0) {
                        mat.ch_rows(i, j);
                        det *= -1;
                        flag = false;
                        break;
                    }
                }
                if (flag) return 0;
            }
            double a = mat.at(i, i);
            for (int j = i+1; j < n; ++j) {
                if (mat.at(j, i) != 0) {
                    double b = mat.at(j, i);
                    for (int k = i; k < m; ++k) mat.at(j, k) -= mat.at(i, k)*(b/a);
                }
            }
        }
        for (int i = 0; i < n; ++i) det *= mat.at(i, i);
        return det;
    }
}

Matrix Matrix::inv() const {
    if (n != m) throw "Matrix::inv - cannot inverse non-square matrix";
    if (!det()) throw "Matrix::inv - degenerate matrix";
    Matrix sys(n, 2*m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            sys.at(i, j) = at(i, j);
        }
        sys.at(i, m+i) = 1;
    }
    for (int i = 0; i < n; ++i) {  // прямой ход
        if (sys.at(i, i) == 0) {
            for (int j = i+1; j < n; ++j) {
                if (sys.at(j, i) != 0) {
                    sys.ch_rows(i, j);
                    break;
                }
            }
        }
        double a = sys.at(i, i);
        for (int j = i+1; j < n; ++j) {
            if (sys.at(j, i) != 0) {
                double b = sys.at(j, i);
                for (int k = i; k < sys.m; ++k) sys.at(j, k) -= sys.at(i, k)*(b/a);
            }
        }
    }
    // обратный ход
    for (int i = n-1; i >= 0; --i) {
        double a = sys.at(i, i);
        for (int k = i; k < sys.m; ++k) sys.at(i, k) /= a;
        for (int j = i-1; j >= 0; --j) {
            double b = sys.at(j, i);
            for (int k = i; k < sys.m; ++k) sys.at(j, k) -= sys.at(i, k)*b;
        }
    }
    Matrix ans(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            ans.at(i, j) = sys.at(i, m+j);
        }
    }
    return ans;
}

Matrix Matrix::to_tr() const {
    Matrix A(n, m, elements);
    int row  = 0;
    for (int col = 0; col < A.m; ++col) {
        if (A.at(row, col) == 0) {
            bool flag = true;
            while (flag && col < A.m) {
                for (int j = row+1; j < n; ++j) {
                    if (A.at(j, col) != 0) {
                        A.ch_rows(row, j);
                        flag = false;
                        break;
                    }
                } if (flag) ++col;
            }
        }
        double a = A.at(row, col);
        for (int j = row+1; j < n; ++j) {
            if (A.at(j, col) != 0) {
                double b = A.at(j, col);
                for (int k = col; k < A.m; ++k) A.at(j, k) -= A.at(row, k)*(b/a);
            }
        }
        ++row;
    }
    return A;
}

int Matrix::rang() const {
    Matrix A(to_tr());
    int i;
    for (i = 0; i < A.n; ++i) {
        bool f = false;
        for (int j = i; j < A.m; ++j) {
            if (A.at(i, j) != 0) {
                f = true;
                break;
            }
        }
        if (!f) break;
    }
    return i;
}

/*Matrix Matrix::dinv() const {
    //if (DetGauss()) return SolveGauss();
    //else {
        Matrix deg_inv(m, n);
        Matrix lambda(m, m);
        for (int i = 1; i < m; ++i) lambda.at(i, i) = 10;
        Matrix A(n ,m, elements);
        deg_inv = ((A.transpose()*A + lambda).inv())*(A.transpose());
        return deg_inv;
    //}
}*/

Matrix EMatrix(int n) {
    Matrix A(n, n);
    for (int i = 0; i < n; ++i) {
        A.at(i, i) = 1;
    }
    return A;
}

Matrix operator+(double i, const Matrix& A) {
    Matrix M(A);
    for (int i = 0; i < M.n; ++i) {
        for (int j = 0; j < M.n; ++j) M.at(i, j) = i+M.at(i, j);
    }
    return M;
}

Matrix operator+(const Matrix& A, double i) {return i+A;}
Matrix operator-(double i, const Matrix& A) {return i+(A*(-1));}
Matrix operator-(const Matrix& A, double i) {return A+(i*(-1));}
Matrix operator*(double i, const Matrix& A) {return A*i;}

Matrix operator/(double i, const Matrix& A) {
    Matrix M(A);
    for (int i = 0; i < M.n; ++i) {
        for (int j = 0; j < M.n; ++j) {
            if (!M.at(i, j)) throw "division by zero";
            M.at(i, j) = 1/M.at(i, j);
        }
    }
    return M;
}

Matrix exp(const Matrix& A) {
    Matrix M(A);
    for (int i = 0; i < M.n; ++i) {
        for (int j = 0; j < M.n; ++j) M.at(i, j) = exp(M.at(i, j));
    }
    return M;
}

Matrix log(const Matrix& A) {
    Matrix M(A);
    for (int i = 0; i < M.n; ++i) {
        for (int j = 0; j < M.n; ++j) M.at(i, j) = log(M.at(i, j));
    }
    return M;
}

double sum(const Matrix& A) {
    double Sum = 0;
    for (int i = 0; i < A.n; ++i) {
        for (int j = 0; j < A.m; ++j) Sum += A.at(i, j);
    }
    return Sum;
}

istream& operator >> (istream& is, Matrix& mat)
{
	for (int row = 0; row < mat.n; ++row)
		for (int col = 0; col < mat.m; ++col)
			is >> mat.at(row, col);
	return is;
}

ostream& operator << (ostream& os, const Matrix& mat)
{
	//os << mat.n << ' ' << mat.m << endl;

	for (int row = 0; row < mat.n; ++row)
	{
		for (int col = 0; col < mat.m; ++col)
			os << mat.at(row, col) << ' ';
		os << endl;
	}
	return os;
}
