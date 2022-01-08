#include <model/finite_difference_approximation.hpp>

#include <cmath>

#include <algorithm/max_diff.hpp>
#include <math/utils.hpp>


namespace {

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
		, h_x_{a_ / static_cast<Float>(N_x_ - 1)}
		, h_z_{b_ / static_cast<Float>(N_z_ - 1)}
	{
	}

	Dependency solve() {
	    Dependency result;

        auto y_k = Matrix(N_x_, std::vector(N_z_, u_0_));
        auto y_k_plus_1 = y_k;
        auto y_k_max = u_0_;
        auto y_k_plus_1_max = 0.0f;

        for (size_t k = 0; /* no condition */; ++k) {
            for (size_t i = 1; i < N_x_ - 1; ++i) {
                const auto xi = static_cast<float>(i) * h_x_;
                for (size_t j = 1; j < N_z_ - 1; ++j) {
                    const auto zj = static_cast<float>(j) * h_z_;
                    y_k_plus_1[i][j] = (y_k[i + 1][j] + y_k[i - 1][j] + y_k[i][j + 1] + y_k[i][j - 1] + f_(xi, zj) / k_ * h_x_ * h_z_) / 4;
                    y_k_plus_1_max = std::max(std::abs(y_k_plus_1[i][j]), y_k_plus_1_max);
                }
            }

            if (std::abs(y_k_plus_1_max - y_k_max) <= eps_) {
                result.t.push_back(static_cast<Float>(k));
                result.y.push_back(y_k_plus_1);
                break;
            }

            std::swap(y_k, y_k_plus_1);

            y_k_max = y_k_plus_1_max;
            y_k_plus_1_max = 0.0f;
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

	const Float h_x_, h_z_;

	Float f_(Float x, Float z) const {
		return f_0_ * std::exp(beta_ * sqr(x - a_ / 2.0f) * sqr(z - b_ / 2.0f));
	}
};

} // namespace


namespace L {

Dependency solve(const Parameters& parameters) {
	return Model(parameters).solve();
}

} // namespace L
