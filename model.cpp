#include "model.h"

#include <cmath>
#include <limits>

#include "integration.h"

#include "quadpackcpp/src/include/workspace.hpp"


ld integrate(function<ld(ld)> fun, ld x1, ld x2) {
	static Workspace<double> workspace(128, 5);

		ld res = 0;

		constexpr double epsabs = 1e-7;
		double abserr;

		workspace.qag(fun, x1, x2, epsabs, 0, res, abserr);

		return res;
}

Model::Model(function<ld(ld, ld)> dzdx_, function<ld(ld, ld)> dzdy_,
			 function<ld(ld)> y0_, function<ld(ld, ld)> beta_, ld alpha_, ld l_)
	: dzdx(dzdx_), dzdy(dzdy_), y0(y0_), beta(beta_), alpha(alpha_), l(l_) {

	int n = ceil(l / dx);

	Points.resize(n + 1);
	sumJ.resize(n);
	sumJ[0] = 0;

	for (int i = 0; i <= n; ++i) {
		Points[i] = {i * l / n, y0(i * l / n)};
	}
}

void Model::optimize(size_t iters) {
	int n = Points.size();


	for (int iter = 0; iter < iters; ++iter) {

		vector<vector<int>> dp_points (2*k + 1, vector<int>(n-1));
		vector<ld> dp_values_prev (2*k + 1, 0);
		vector<ld> dp_values_cur (2*k + 1, 0);
		vector<ld> add_prev(2*k+1, 0);
		vector<ld> add_cur(2*k+1, 0);

		for (int t = 0; t < n-1; ++t) {

			auto [x1,y1] = Points[t];
			auto [x2,y2] = Points[t+1];

			for (int i = -k; i <= k; ++i) {
				if (t == n-2 && i != 0) continue;

				ld J_local_best = numeric_limits<ld>::max();
				int best_k = 0;
				ld best_a;
				ld best_b;

				for (int j = -k; j <= k; ++j) {
					if (t == 0 && j != 0) continue;
					
					y1 += dy * j;
					y2 += dy * i;

					ld a = (y2-y1) / dx;
					ld b = y1 - a*x1;

					ld J_local_cur = J(x1, x2, a, b, add_prev[j+k]) + dp_values_prev[j+k];

					if (J_local_cur < J_local_best) {
						J_local_best = J_local_cur;
						best_k = j;
						best_a = a;
						best_b = b;
					}
				}

				add_cur[i+k] = add_prev[best_k + k] + helper(x1, x2, best_a, best_b);
				dp_values_cur[i+k] = J_local_best;
				dp_points[i+k][t] = best_k + k; 
			}

			swap(dp_values_prev, dp_values_cur);
			swap(add_prev, add_cur);
		}

		J_total_value = dp_values_prev[k];
	
		vector<Point> new_Points;
		new_Points.resize(n);

		new_Points[n-1] = Points[n-1]; 
		int cur = k;

		for (int t = n-2; t >= 0; --t)
		{
			int j = dp_points[cur][t];
			auto [x,y] = Points[t];
			y += (j-k) * dy;
			new_Points[t] = {x,y};
			cur = j;
		}

		Points = move(new_Points);
	}
}
/*
ld Model::J() {
	int n = Points.size();

	ld out = 0;

	for (int i = 0; i < n - 1; ++i) {
		ld a =
			(Points[i + 1].y - Points[i].y) / (Points[i + 1].x - Points[i].x);
		ld b = Points[i].y - a * Points[i].x;

		out += J(Points[i].x, Points[i + 1].x, a, b);
	}

	return out;
}*/

ld Model::helper(ld x1, ld x2, ld a, ld b) {
	function<ld(ld)> f1{
		[a, b, dzdx = this->dzdx, dzdy = this->dzdy](ld x) -> ld {
			ld y = a * x + b;
			ld dzdx_value = dzdx(x, y) + dzdy(x, y) * a;

			return sqrt(1 + a*a + pow(dzdx_value, 2));
	}};

	return integrate(f1,x1,x2);
}
 
 
ld Model::J(ld x1, ld x2, ld a, ld b, ld add) {
	function<ld(ld)> f1{
		[a, b, dzdx = this->dzdx, dzdy = this->dzdy](ld x) -> ld {
			ld y = a * x + b;
			ld dzdx_value = dzdx(x, y) + dzdy(x, y) * a;

			return sqrt(1 + a*a + pow(dzdx_value, 2));
	}};

	function<ld(ld)> f2{
		[a, b, dzdx = this->dzdx, dzdy = this->dzdy, beta = this->beta](ld x) -> ld {
			ld y = a * x + b;
			ld dzdx_value = dzdx(x, y) + dzdy(x, y) * a;

			return beta(x,y) * sqrt(1 + a*a + pow(dzdx_value, 2));
	}};

	ld integr1 = integrate(f1, x1, x2);
	ld res1 = (add + 0.5 * integr1) * integr1 * alpha; 
	ld res2 = integrate(f2, x1, x2);

	return res1 + res2;
}

ld Model::dx{0.01};
ld Model::dy{0.01};
int Model::k{4};