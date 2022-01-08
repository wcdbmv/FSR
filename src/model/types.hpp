#pragma once

#include <math/types.hpp>


struct Parameters {
	Float a, b;

	Size N_x, N_z;
	Float tau;
	Float eps;

	Float u_0;

	Float f_0;
	Float beta;

	Float k;

	Size N_p = 10;
};

struct Dependency {
	Vector t;
	Container<Matrix> y;
};
