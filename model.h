#pragma once

#include <functional>
#include <map>
#include <vector>

using std::function, std::map, std::vector;
using ld = double;

ld integrate(function<ld(ld)> fun, ld x1, ld x2);

struct Point {
	ld x;
	ld y;
};

struct Model {
	Model(function<ld(ld, ld)> dzdx_, function<ld(ld, ld)> dzdy_,
		  function<ld(ld)> y0, function<ld(ld, ld)> beta_, ld alpha_, ld l_);

	void optimize(size_t iter);

	// y = ax + b
	//ld J();
	ld J(ld x1, ld x2, ld a, ld b, ld add);
	ld helper (ld x1, ld x2, ld a, ld b);

	// calculated y:
	ld J_total_value = -1;
	vector<Point> Points;
	vector<ld> sumJ;

	// task parameters
	function<ld(ld, ld)> dzdx;
	function<ld(ld, ld)> dzdy;
	function<ld(ld, ld)> beta;
	function<ld(ld)> y0;

	ld alpha;
	ld l;

	// method parameters
	static ld dx;
	static ld dy;
	static int k;
};
