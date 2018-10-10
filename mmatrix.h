#ifndef MMATRIX_H_INCLUDED
#define MMATRIX_H_INCLUDED

#include<iostream>
#include<cmath>

class Matrix {
protected:
    double* elements;
    void ch_rows(int row1, int row2);
public:
	const int n;
	const int m;

	Matrix() : elements(new double[1]), n(1), m(1) {elements[0] = 0;}
	Matrix(const Matrix&);
	Matrix(int, int);
	Matrix(const int, const int, double*);

	double& at(int, int);
	double at(int, int) const;

	void operator=(const Matrix&);
	Matrix transpose() const;

	Matrix operator*(double) const;
	Matrix operator/(double) const;
	Matrix operator+(const Matrix&) const;
	Matrix operator-(const Matrix&) const;
	Matrix operator*(const Matrix&) const;
	double operator()(int i, int j) {return at(i, j);}

	Matrix to_tr() const;
	int rang() const;
	double det() const;
	Matrix inv() const;
	//Matrix dinv() const;

	~Matrix() { delete[] elements; }

};

Matrix EMatrix(int);

Matrix operator+(double, const Matrix&);
Matrix operator+(const Matrix&, double);
Matrix operator-(double, const Matrix&);
Matrix operator-(const Matrix&, double);
Matrix operator*(double, const Matrix&);
Matrix operator/(double, const Matrix&);
Matrix log(const Matrix&);
Matrix exp(const Matrix&);
double sum(const Matrix&);

std::istream& operator >> (std::istream& is, Matrix& mat);
std::ostream& operator << (std::ostream& os, const Matrix& mat);

#endif // MMATRIX_H_INCLUDED
