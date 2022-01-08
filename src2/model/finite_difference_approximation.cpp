#include <model/finite_difference_approximation.hpp>

#include <cmath>

#include <algorithm/max_diff.hpp>
#include <algorithm/thomas_algorithm.hpp>
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
		, h_x_2_{sqr(h_x_)}
		, h_z_2_{sqr(h_z_)}
	{
	}

	Dependency solve() {
		auto y_prev = Matrix(N_x_, Vector(N_z_, u_0_));
		auto y_mid = y_prev;
		auto y_curr = y_prev;

		Dependency result;
		result.t.push_back(0.0f);
		result.y.push_back(y_curr);

		for (auto t = tau_; t <= 700.0f; t += tau_) {
			for (Size j = 1; j < N_z_ - 1; ++j) {
				const auto A = A_for_overline_y_();
				const auto B = B_for_overline_y_();
				const auto C = C_for_overline_y_();
				const auto D = D_for_overline_y_(j, y_prev);

				const auto y_ij_mid = ThomasAlgorithm(A, B, C, D);

				for (Size i = 0; i < N_x_; ++i) {
					y_mid[i][j] = y_ij_mid[i];
				}
			}

			for (Size i = 1; i < N_x_ - 1; ++i) {
				const auto A = A_for_hat_y_();
				const auto B = B_for_hat_y_();
				const auto C = C_for_hat_y_();
				const auto D = D_for_hat_y_(i, y_mid);

				ThomasAlgorithm(A, B, C, D, y_curr[i]);
			}

			result.t.push_back(t);
			result.y.push_back(y_curr);

			if (max_diff(y_prev, y_curr) < eps_) {
				break;
			}

			auto tmp = std::move(y_prev);
			y_prev = std::move(y_curr);
			y_curr = std::move(tmp);
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
	const Float h_x_2_, h_z_2_;

	Float f_(Float x, Float z) const {
		return f_0_ * std::exp(beta_ * sqr(x - a_ / 2.0f) * sqr(z - b_ / 2.0f));
	}

	const Float A_i_for_overline_y_ = tau_ * h_z_2_;
	const Float B_i_for_overline_y_ = -2 * h_z_2_ * (tau_ + h_x_2_);
	const Float C_i_for_overline_y_ = tau_ * h_z_2_;

	Float D_i_for_overline_y_(Float xi, Float zj, Float y_i_jm1, Float y_i_j, Float y_i_jp1) const {
		return h_x_2_ * (2 * (tau_ - h_z_2_) * y_i_j - tau_ * (y_i_jm1 + y_i_jp1 + h_z_2_ * f_(xi, zj) / k_));
	}

	Vector A_for_overline_y_() const {
		Vector A(N_x_, A_i_for_overline_y_);
		A.front() = 0;
		A.back() = 0;
		return A;
	}

	Vector B_for_overline_y_() const {
		Vector B(N_x_, B_i_for_overline_y_);
		B.front() = 1;
		B.back() = 1;
		return B;
	}

	Vector C_for_overline_y_() const {
		Vector C(N_x_, C_i_for_overline_y_);
		C.front() = 0;
		C.back() = 0;
		return C;
	}

	Vector D_for_overline_y_(Size j, const Matrix& y_prev) const {
		const auto zj = static_cast<Float>(j) * h_z_;

		Vector D(N_x_);
		D.front() = u_0_;
		for (Size i = 1; i < N_x_ - 1; ++i) {
			const auto xi = static_cast<Float>(i) * h_x_;
			D[i] = D_i_for_overline_y_(xi, zj, y_prev[i][j - 1], y_prev[i][j], y_prev[i][j + 1]);
		}
		D.back() = u_0_;

		return D;
	}

	const Float A_j_for_hat_y_ = tau_ * h_x_2_;
	const Float B_j_for_hat_y_ = -2 * h_x_2_ * (tau_ + h_z_2_);
	const Float C_j_for_hat_y_ = tau_ * h_x_2_;

	Float D_j_for_hat_y_(Float xi, Float zj, Float y_im1_j, Float y_i_j, Float y_ip1_j) const {
		return h_z_2_ * (2 * (tau_ - h_x_2_) * y_i_j - tau_ * (y_im1_j + y_ip1_j + h_x_2_ * f_(xi, zj) / k_));
	}

	Vector A_for_hat_y_() const {
		Vector A(N_z_, A_j_for_hat_y_);
		A.front() = 0;
		A.back() = 0;
		return A;
	}

	Vector B_for_hat_y_() const {
		Vector B(N_z_, B_j_for_hat_y_);
		B.front() = 1;
		B.back() = 1;
		return B;
	}

	Vector C_for_hat_y_() const {
		Vector C(N_z_, C_j_for_hat_y_);
		C.front() = 0;
		C.back() = 0;
		return C;
	}

	Vector D_for_hat_y_(Size i, const Matrix& y_mid) const {
		const auto xi = static_cast<Float>(i) * h_x_;

		Vector D(N_z_);
		D.front() = u_0_;
		for (Size j = 1; j < N_z_ - 1; ++j) {
			const auto zj = static_cast<Float>(j) * h_z_;
			D[j] = D_j_for_hat_y_(xi, zj, y_mid[i - 1][j], y_mid[i][j], y_mid[i + 1][j]);
		}
		D.back() = u_0_;

		return D;
	}
};

} // namespace


namespace FDA {

Dependency solve(const Parameters& parameters) {
	return Model(parameters).solve();
}

} // namespace FDA
