#include <cmath>
#include <compare>
#include <fstream>
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
		[](ld x) -> ld { return  x - 3.775e-1 * sin(3.14*x) + 6.80e-2*sin(3.14*2*x) - 7.967e-4 * sin(3.14*3*x); }};

	function<ld(ld, ld)> beta{[](ld x, ld y) -> ld { return 0.5; }};

	Model model(dzdx, dzdy, y0, beta, alpha, l);

	model.optimize(10);


	std::cerr << model.J_total_value << '\n';

	for (auto [x, y] : model.Points) {
		cout << x << ' ' << y << '\n';
	}
}