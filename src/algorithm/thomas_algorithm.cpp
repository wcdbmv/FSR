#include <algorithm/thomas_algorithm.hpp>

#include <stdexcept>
#include <tuple>


namespace {

std::tuple<Vector, Vector> ForwardSweep(const Vector& a, const Vector& b, const Vector& c, const Vector& d) {
	const auto N = a.size();

	Vector c_star(N);
	Vector d_star(N);

	c_star[0] = c[0] / b[0];
	d_star[0] = d[0] / b[0];

	for (Size i = 1; i < N; ++i) {
		const auto denominator = b[i] - a[i] * c_star[i - 1];
		c_star[i] = c[i] / denominator;
		d_star[i] = (d[i] - a[i] * d_star[i - 1]) / denominator;
	}

	return {c_star, d_star};
}

void ReverseSweep(const Vector& c_star, const Vector& d_star, Vector& x) {
	const auto N = c_star.size();

	x[N - 1] = d_star[N - 1];

	for (int i = static_cast<int>(N) - 2; i >= 0; --i) {
		const auto j = static_cast<Size>(i);
		x[j] = d_star[j] - c_star[j] * x[j + 1];
	}
}

} // namespace


Vector ThomasAlgorithm(const Vector& a, const Vector& b, const Vector& c, const Vector& d) {
	Vector x(a.size());
	ThomasAlgorithm(a, b, c, d, x);
	return x;
}

void ThomasAlgorithm(const Vector& a, const Vector& b, const Vector& c, const Vector& d, Vector& x) {
	const auto N = a.size();

	if (b.size() != N || c.size() != N || d.size() != N || x.size() != N || !N) {
		throw std::invalid_argument("ThomasAlgorithm: vectors must be the same size and non-empty");
	}

	if (a.front() != 0 || c.back() != 0) {
		throw std::invalid_argument("ThomasAlgorithm: a.front() and c.back() must be zero");
	}

	const auto [c_star, d_star] = ForwardSweep(a, b, c, d);

	ReverseSweep(c_star, d_star, x);
}
