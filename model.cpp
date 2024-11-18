#include "model.h"

#include <cmath>
#include <limits>

#include "integration.h"

ld integrate(function<ld(ld)> fun, ld x1, ld x2) {
	ld result = biv::makeCompoundSF(x1, x2, fun);
	return result;
}

Model::Model(function<ld(ld, ld)> dzdx_, function<ld(ld, ld)> dzdy_,
			 function<ld(ld)> y0_, function<ld(ld, ld)> beta_, ld alpha_, ld l_)
	: dzdx(dzdx_), dzdy(dzdy_), y0(y0_), beta(beta_), alpha(alpha_), l(l_) {
	int n = ceil(l / dx);

	Points.resize(n + 1);

	for (int i = 0; i <= n; ++i) {
		Points[i] = {i * l / n, y0(i * l / n)};
	}

	for (int i = 0; i < n; ++i) {
		auto [x1,y1] = Points[i];
		auto [x2,y2] = Points[i+1];

		ld a = (y2-y1) / dx;
		ld b = y1 - a*x1;

		J_total_value += J(x1, x2, a, b);
	}
}

void Model::optimize(size_t iters) {
	int n = Points.size();

	for (int iter = 0; iter < iters; ++iter) {

		vector<vector<int>> dp_points (2*k + 1, vector<int>(n-1));
		vector<ld> dp_values_prev (2*k + 1, 0);
		vector<ld> dp_values_cur (2*k + 1, 0);

		for (int t = 0; t < n-1; ++t) {

			for (int i = -k; i <= k; ++i) {
				if (t == n-2 && i != 0) continue;

				ld J_local_best = numeric_limits<ld>::max();
				int best_k = 0;

				for (int j = -k; j <= k; ++j) {
					if (t == 0 && j != 0) continue;

					auto [x1,y1] = Points[t];
					auto [x2,y2] = Points[t+1];
					
					y1 += dy * j;
					y2 += dy * i;

					ld a = (y2-y1) / dx;
					ld b = y1 - a*x1;

					ld J_local_cur = J(x1, x2, a, b) + dp_values_prev[j+k];

					if (J_local_cur < J_local_best) {
						J_local_best = J_local_cur;
						best_k = j;
					}
				}

				dp_values_cur[i+k] = J_local_best;
				dp_points[i+k][t] = best_k + k; 
			}

			swap(dp_values_prev, dp_values_cur);
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

ld Model::J() {
	int n = Points.size() - 1;

	ld out = 0;

	for (int i = 0; i < n; ++i) {
		ld a =
			(Points[i + 1].y - Points[i].y) / (Points[i + 1].x - Points[i].x);
		ld b = Points[i].y - a * Points[i].x;

		out += J(Points[i].x, Points[i + 1].x, a, b);
	}

	return out;
}

ld Model::J(ld x1, ld x2, ld a, ld b) {
	function<ld(ld)> f1{
		[a, b, dzdx = this->dzdx, dzdy = this->dzdy](ld x) -> ld {
			ld y = a * x + b;
			ld dzdx_value = dzdx(x, y) + dzdy(x, y) * a;

			return sqrt(1 + pow(a, 2) + pow(dzdx_value, 2));
		}};

	function<ld(ld)> f2{[a, b, dzdx = this->dzdx, dzdy = this->dzdy,
						 beta = this->beta](ld x) -> ld {
		ld y = a * x + b;
		ld dzdx_value = dzdx(x, y) + dzdy(x, y) * a;

		return beta(x, y) * sqrt(1 + pow(a, 2) + pow(dzdx_value, 2));
	}};

	ld res =
		alpha * 0.5 * pow(integrate(f1, x1, x2), 2) + integrate(f2, x1, x2);
	return res;
}

ld Model::dx{0.02};
ld Model::dy{0.02};
int Model::k{2};
