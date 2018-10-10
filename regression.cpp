#include "regretion.h"

using namespace std;

double costFunction(const Matrix& X, const Matrix& y, const Matrix& theta, double lambda) {
	return ((((X*theta - y).transpose())*(X*theta - y))*((double)1/(2*(X.n)))+(theta.transpose()*theta)*(lambda/(2*X.n))).at(0, 0);
}

Matrix& gradient(const Matrix& X, const Matrix& y, Matrix& theta, double alpha, double lambda, int iterations_number) {
	for (int i = 0; i < iterations_number; ++i) {
        //cout << costFunction(X, y, theta, lambda) << endl;
        theta = theta - (((X.transpose())*(X*theta - y))*(alpha/(double)X.n)) - (theta*lambda);
	}
	return theta;
}

Matrix& normalEquasion(const Matrix& X, const Matrix& y, Matrix& theta, double lambda) {
    theta = (((X.transpose())*X).inv())*(X.transpose())*y;
    return theta;
}

void Normalizator::clear() {
    delete[] MAX;
    delete[] MIN;
    delete[] SUM;
    MAX = NULL;
    MIN = NULL;
    SUM = NULL;
    size = 0;
}

double Normalizator::normFeature(double x, int feature_num) {return (x - SUM[feature_num])/(MAX[feature_num] - MIN[feature_num]);}

Matrix Normalizator::normVector(const Matrix& X) {
    Matrix X_norm(X);
    for (int i = 1; i < size; ++i) X_norm.at(0, i) = normFeature(X_norm.at(0, i), i);
    return X_norm;
}

Matrix Normalizator::operator()(const Matrix& X) {
    size = X.m;
    MIN = new double[X.m];
    MAX = new double[X.m];
    SUM = new double[X.m];
    Matrix XN(X);
	for (int i = 1; i < XN.m; ++i) {
		for (int j = 0; j < XN.n; ++j) {
			if (j == 0) {
				MAX[i] = XN.at(0, i);
				MIN[i] = XN.at(0, i);
				SUM[i] = XN.at(0, i);
			} else {
				if (XN.at(j, i) > MAX[i]) MAX[i] = XN.at(j, i);
				if (XN.at(j, i) < MIN[i]) MIN[i] = XN.at(j, i);
				SUM[i] += XN.at(j, i);
			}
		}
		SUM[i] /= XN.n;
		for (int j = 0; j < XN.n; ++j) {
			XN.at(j, i) = normFeature(XN.at(j, i), i);
		}
	}
	return XN;
}
