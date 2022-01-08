#include <model/finite_difference_approximation.hpp>

#include <cmath>
#include <random>

#include <algorithm/max_diff.hpp>
#include <math/utils.hpp>


namespace {

template <int Min, int Max>
int randint() {
	static thread_local std::mt19937 generator(std::random_device{}());
	static thread_local std::uniform_int_distribution<int> distribution(Min, Max);
	return distribution(generator);
}

enum class Direction {
	Left = 0,
	Right = 1,
	Up = 2,
	Down = 3,
};

Direction random_direction() {
	return static_cast<Direction>(randint<0, 3>());
}

Float mean(const Vector& vector) {
	return std::accumulate(vector.begin(), vector.end(), Float{}) / static_cast<Float>(vector.size());
}

class Model {
public:
	explicit Model(const Parameters& parameters)
		: a_{parameters.a}
		, b_{parameters.b}
		, N_x_{parameters.N_x}
		, N_z_{parameters.N_z}
		, tau_{parameters.tau}
		, eps_{parameters.eps}
		, u_0_{parameters.u_0}
		, f_0_{parameters.f_0}
		, beta_{parameters.beta}
		, k_{parameters.k}
		, N_p_{parameters.N_p}
		, h_x_{a_ / static_cast<Float>(N_x_ - 1)}
		, h_z_{b_ / static_cast<Float>(N_z_ - 1)}
	{
	}

	Dependency solve() {
		Dependency result;

		Matrix y(N_x_, Vector(N_z_, u_0_));

		const auto final_layer = std::min((N_x_ - 1) / 2, (N_z_ - 1) / 2);

		for (Size layer = 1; layer <= final_layer; ++layer) {
			for (Size j = layer; j + layer < N_z_; ++j) {
				y[layer][j] = simulate_n(layer, j, layer - 1, y);
				y[N_x_ - layer - 1][j] = simulate_n(N_x_ - layer - 1, j, layer - 1, y);
			}
			for (Size i = layer + 1; i + layer + 1 < N_x_; ++i) {
				y[i][layer] = simulate_n(i, layer, layer - 1, y);
				y[i][N_z_ - layer - 1] = simulate_n(i, N_z_ - layer - 1, layer - 1, y);
			}
			result.t.push_back(static_cast<Float>(layer));
			result.y.push_back(y);
		}

		return result;
	}

private:
	const Float a_, b_;

	const Size N_x_, N_z_;
	const Float tau_;
	const Float eps_;

	const Float u_0_;

	const Float f_0_;
	const Float beta_;

	const Float k_;

	const Size N_p_;

	const Float h_x_, h_z_;

	Float x_i_(Size i) const {
		return static_cast<Float>(i) * h_x_;
	}

	Float z_j_(Size j) const {
		return static_cast<Float>(j) * h_z_;
	}

	Float f_(Float x, Float z) const {
		return f_0_ * std::exp(beta_ * sqr(x - a_ / 2.0f) * sqr(z - b_ / 2.0f));
	}

	Float f_(Size i, Size j) const {
		return f_(x_i_(i), z_j_(j));
	}

	Float simulate(Size i, Size j, Size border, const Matrix& y) {
		Vector fs{f_(i, j)};

		Size steps = 0;
		while (border < i && i < N_x_ - border - 1 && border < j && j < N_z_ - border - 1) {
			const auto direction = random_direction();
			switch (direction) {
			case Direction::Left:
				--j;
				break;
			case Direction::Right:
				++j;
				break;
			case Direction::Up:
				--i;
				break;
			case Direction::Down:
				++i;
				break;
			}
			fs.push_back(f_(i, j));
			++steps;
		}

		return y[i][j] + mean(fs) * h_x_ * h_z_ / (4 * k_) * static_cast<Float>(steps);
	}

	Float simulate_n(Size i, Size j, Size border, const Matrix& y) {
		const auto N = N_p_ / sqr(border + 1);

		Vector ys;
		ys.reserve(N);

		for (Size k = 0; k < N; ++k) {
			ys.push_back(simulate(i, j, border, y));
		}

		return mean(ys);
	}
};

} // namespace


namespace P {

Dependency solve(const Parameters& parameters) {
	return Model(parameters).solve();
}

} // namespace P
