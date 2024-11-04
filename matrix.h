#pragma once
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace std;
namespace biv {
template <typename Number>
class Matrix;
template <typename Number>
ostream &operator<<(ostream &os, const Matrix<Number> &a);

template <typename Number>
class Matrix {
   public:
	vector<vector<Number>> v;
	size_t n;
	size_t m;

	Matrix() : n(0), m(0) {}
	Matrix(const vector<vector<Number>> &v)
		: v(v), n(v.size()), m(v[0].size()) {}
	Matrix(const vector<Number> &v) : n(1), m(v.size()) {
		vector<vector<Number>> copy(1, vector<Number>(v.size()));
		copy[0] = v;
		this->v = copy;
	}
	Matrix(const initializer_list<initializer_list<Number>> &init) {
		n = init.size();
		m = (*init.begin()).size();
		for (auto data : init) {
			v.push_back(data);
		}
	}
	Matrix(const initializer_list<Number> &init)
		: v{init}, n(1), m(init.size()) {}
	Matrix(size_t n, size_t m) : n(n), m(m) {
		vector<vector<Number>> q(n, vector<Number>(m));
		v = q;
	}
	Matrix(size_t n) : n(n), m(n) {	 // make E
		vector<vector<Number>> q(n, vector<Number>(n, 0));
		for (int i = 0; i < n; i++) q[i][i] = 1;
		v = q;
	}
	Matrix(size_t n, size_t m, Number t) : n(n), m(m) {
		vector<vector<Number>> q(n, vector<Number>(m, t));
		v = q;
	}
	// copy constr
	Matrix(const Matrix &a) {
		n = a.n;
		m = a.m;
		v = a.v;
	}
	Number &operator()(int i, int j) {
		if (i >= this->n || j >= this->m || i < 0 || j < 0) {
			cout << "list of matrix is out of range\n";
			throw;
		}
		return this->v[i][j];
	}

	const Number &operator()(int i, int j) const {
		if (i >= this->n || j >= this->m || i < 0 || j < 0) {
			cout << "Seg Fault\n";
			throw;
		}
		return this->v[i][j];
	}
	vector<Number> &operator[](int i) {
		if (i >= this->n || i < 0) {
			cout << "list of matrix is out of range\n";
			throw;
		}
		return this->v[i];
	}
	// Multiply by number
	Matrix operator*(Number x) {
		Matrix B = *this;
		for (auto &i : B.v)
			for (auto &j : i) j *= x;
		return B;
	}
	Matrix operator+(Matrix A) {
		Matrix<Number> res = *this;
		if (this->n != A.n || this->m != A.m) {
			cout << "bad sum of matrices\n";
			throw;
		}
		for (int i = 0; i < A.n; i++) {
			for (int j = 0; j < A.m; j++) res(i, j) += A(i, j);
		}
		return res;
	}

	Matrix &operator*=(Number x) {
		Matrix m = this->v;
		for (auto &i : this->v)
			for (auto &j : i) j *= x;
		return m;
	}

	Matrix &operator+=(Matrix a) {
		*this = *this + a;
		return *this;
	}
	Matrix &operator=(Matrix a) {
		n = a.n;
		m = a.m;
		v = a.v;
		return *this;
	}

	Matrix operator^(int k) {
		if (this->n != this->m) {
			cout << "n != m \n";
			throw;
		}
		Matrix<Number> Res(this->n, this->n, static_cast<Number>(0));
		Matrix<Number> A = *this;
		for (int i = 0; i < n; i++) Res(i, i) = 1;
		while (k > 0) {
			if ((k & 1) != 0) Res = Res * A;
			A = A * A;
			k >>= 1;
		}
		return Res;
	}

	// Matrix Multipl
	Matrix operator*(Matrix b) {
		// cout << this->m << " " << b.n << "<------\n";
		if (this->m != b.n) {
			cout << "Bad dimensions in multiplying\n";
			throw;
		}
		double curr_sum = 0;
		Matrix res(n, b.m);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < b.m; j++) {	 // i, j position in new matrix
				res(i, j) = 0;
				for (int q = 0; q < m; q++) {
					res(i, j) += (*this)(i, q) * b(q, j);
				}
			}
		}
		return res;
	}
	// ÏÎÑÌÎÒÐÅÒÜ ËÅÊÖÈÞ
	/*Matrix& operator*=(Matrix a) {
			*this = *this * a;
			return *this;
	}
  */ Matrix transpose() {
		Matrix &a = *this;
		Matrix b(a.m, a.n);
		for (int i = 0; i < a.m; i++)
			for (int j = 0; j < a.n; j++) b.v[i][j] = a.v[j][i];
		return b;
	}
	Matrix inverse() {
		Matrix<Number> &A = *this;
		if (A.n != A.m) {
			cout << "Matrix hasnt inverse\n";
			throw;
		}
		Matrix<Number> ANS(A.n, A.m, static_cast<Number>(0));

		for (int i = 0; i < A.n; i++) {
			Matrix<Number> B(A.n, 1, static_cast<Number>(0));
			B(i, 0) = 1;
			Matrix<Number> ans = GaussSlau(A, B);
			for (int j = 0; j < A.n; j++) {
				ANS(j, i) = ans(j, 0);
			}
		}
		return ANS;
	}
};
template <typename Number>
ostream &operator<<(ostream &os, const Matrix<Number> &a) {
	for (auto i : a.v) {
		for (double j : i) os << setprecision(10) << j << " ";
		os << '\n';
	}
	return os;
}
template <typename Number>
Matrix<Number> GaussSlau(const Matrix<Number> &A, const Matrix<Number> &b) {
	size_t n = b.n;
	if (A.m != A.n) {
		cout << "Matrix A is not a square!\n";
		throw;
	}
	Matrix<Number> W(n, n + 1);
	for (int i = 0; i < n; i++) {
		copy((A.v[i]).begin(), (A.v[i]).end(), (W.v[i]).begin());
		W.v[i][n] = b.v[i][0];
	}
	// cout << "W:\n" << W;
	for (int j = 0; j < n; j++) {
		int maxx = j;
		for (int i = j + 1; i < n; i++)
			if (abs(W(maxx, j)) < abs(W(i, j))) {
				swap(W.v[maxx], W.v[i]);
				maxx = i;
			}
		double denom = W(j, j);
		for (int i = 0; i < n + 1; i++) {
			W(j, i) = W(j, i) * 1.0 / denom;
		}
		for (int i = j + 1; i < n; i++) {
			denom = W(i, j);
			for (int q = 0; q < n + 1; q++) {
				W(i, q) = W(i, q) - (W(j, q) * denom);
			}
		}
	}
	// cout << "W:\n" << W;
	//  inverse
	for (int j = static_cast<int>(n) - 1; j > 0; j--) {
		// for (int dd = 0; dd < n; dd++)
		//	W.v[j][dd] /= W.v[j][j];

		for (int i = j - 1; i >= 0; i--) {
			double denom = W(i, j);
			for (int dd = 0; dd < n + 1; dd++)	// CAN GET MORE EFFIENCY
				W(i, dd) = W(i, dd) - (W(j, dd) * denom);
		}
	}
	Matrix<Number> res(n, 1);
	for (size_t i = 0; i < n; i++) res.v[i][0] = W.v[i][n];
	return res;
}
}  // namespace biv