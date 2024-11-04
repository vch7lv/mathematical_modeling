#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <numeric>
#include <fstream>
#include <functional>
#include "matrix.h"

// namespace var_b {
// 	static const double alpha = 1.0 / 3;
// 	static const double a = 1.5;
// 	static const double b = 3.3;
// }
/*
static const double alpha = 1;
static const double a = 0;
static const double b = 3;*/

using vpairdd = std::vector<std::pair<long double, long double>>;

vpairdd evenNet(int p, long double x1, long double x2) {
	vector<pair<long double,long double>> res(p);
	double h = (x2 - x1) / p;
	double curr = x1;
	for (int i = 0; i < p; i++) {
		res[i] = make_pair(curr, curr + h);
		curr += h;
	}
	return res;
}

// double func(double x) {
// 	return (2.0 * cos(2.5 * x) * exp(x / 3.0) + 4 * sin(3.5 * x) * exp(-3 * x) + x);
// }
// double func_test(double x) {
// 	return x * x;

namespace biv {
	double computeMoment(double k, double a, double b, long double x1, long double x2);
	template <typename Number>
	double computeIntegral(Matrix<Number> Mu, const vector<double>& nodes);
	//template <typename Number>
	//Matrix<Number> RhimannIntegrale(double, double, int, Matrix<Number>(*MatrixFunc)(double ksi));

	vector<double> makeNodes(int n, double a, double b) {
		vector<double> nodes(n);
		double h = (b - a) / (n-1);
		nodes[0] = a;
		for (int i = 1; i < n; i++)
			nodes[i] = nodes[i - 1] + h;
		return nodes;
	}
	
	template <typename Number>
	Matrix<Number> makeISF(int n, vector<double>& nodes, double a, double b) {
		nodes = makeNodes(n, a, b);
		//for (auto i : nodes) cout << i << ' ';
		//cout << '\n';
		Matrix<Number> Mu(nodes.size(), 1);
		for (int i = 0; i < nodes.size(); i++)
			Mu(i, 0) = computeMoment(i, a, b, );
		Matrix<Number> A(nodes.size(), nodes.size());
		for (int i = 0; i < nodes.size(); i++)
			for (int j = 0; j < nodes.size(); j++)
				A(i, j) = pow(nodes[j] - a, i); //t = x-a   f(xj - a)
		Matrix<Number> Ashki = GaussSlau(A, Mu);
		bool flag = false;
		double Ashki_sum = 0, XD_A_sum = 0;
		for (int i = 0; i < n; i++) {
      		Ashki_sum += abs(Ashki(i, 0));
			XD_A_sum += Ashki(i, 0);
			if (Ashki(i, 0) < 0) flag = true;
		}
		//cout << Ashki_sum << '\n';
		//fout << Ashki_sum << ',' << n << '|';
		//if (flag) cout << "\nthere is some negative A[i]!!!\n";
		return Ashki;
	}
	double computeMoment(double k, double a, double b, long double x1, long double x2) {
		double ans;
		if (k == -1) {
			ans = 1e9;
		}else {
			ans = pow((b - x1), k + 1);
			ans /= k + 1;
			ans -= pow((a - x1), k + 1)/(k+1);
		}
		return ans;
	}
	//important to note rhat nodes must lay between a...b
	template <typename Number>
	double computeIntegral(Matrix<Number> Aj, const vector<double>& nodes,
                           const std::function<long double(long double)>& func) {
		int n = Aj.n;
		double integrale = 0;
		vector<double> sl(n);
		for (int i = 0; i < n; i++) {
			sl[i] = Aj(i, 0) * func(nodes[i]);// xj = tj + a(a a + b / 2 b)
		}
		sort(sl.begin(), sl.end());
		for (auto i : sl) {
			integrale += i;
		}
		return integrale;
	}

	// n - ?????????�???�???? ??????-???? ??????????, p - ??????-???? ???�?????????�??????, rule - ???????????? ??????????-?????�????????
	long double makeCompoundSF(long double x1, long double x2,
                               const std::function<long double(long double)>& func) {
        int n = 50;
        int p = 5;
        function<vpairdd(int, long double, long double)> rule = evenNet;
		vpairdd interv = rule(p, x1, x2);
		//???????????? ?�?????????�-????????
		vector<vector<double>> nodes(p);
		vector<Matrix<double>> Ashki(p);
		vector<int> node_counter(p, 0);
		if (interv.size() != p) {
			throw runtime_error("");
		}
		//???????????�?�???? ?�????????, ???????????? ???� ?????�?????�?? ???????? ???????� ?�?�???? ??????????????
		while (n > 0) {
			for (int i = 0; i < p; i++) {
				if (n > 0) {
					node_counter[i]++;
					n--;
				}
				else break;
			}
		}
		for (int i = 0; i < p; i++) {
			if (i != p - 1) {
				Ashki[i] = makeISF<double>(node_counter[i], nodes[i], interv[i].first, interv[i].second);
			}
			else {
				Ashki[i] = makeISF<double>(node_counter[i], nodes[i], interv[i].first, interv[i].second);
			}
		}
		//make sure that we do it right
		double overall_curr = 0;
		for (int i = 0; i < p; i++) {
			double curr = 0;
			for (int j = 0; j < node_counter[i]; j++) {
				curr += Ashki[i](j, 0);
			}
			overall_curr += curr;
		}
		double ans = 0;
		for (int i = 0; i < p; i++) {
			double XD = computeIntegral<double>(Ashki[i], nodes[i], func);
			ans += XD;
		}
		return ans;
	}

}