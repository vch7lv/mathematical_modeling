#include <cmath>
#include <compare>
#include <iostream>

#include "model.h"

using std::cout;

int main() {
	ld alpha = 0.1;
	ld l = 1;

	function<ld(ld, ld)> dzdx{
		[](ld x, ld y) -> ld { return 5 * sin(y) * cos(5 * x); }};

	function<ld(ld, ld)> dzdy{
		[](ld x, ld y) -> ld { return sin(5 * x) * cos(y); }};

	function<ld(ld)> y0{
		[](ld x) -> ld { return x; }};

	function<ld(ld, ld)> beta{[](ld x, ld y) -> ld { return 0.5; }};

	Model model(dzdx, dzdy, y0, beta, alpha, l);

	//model.optimize(20);

	cout << model.J(0,1,1,0) << '\n';
	//cout << model.J_total_value << '\n';

	for (auto [x, y] : model.Points) {
		cout << x << ' ' << y << '\n';
	}
}