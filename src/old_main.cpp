#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// Вариант 3. Математическая модель с постоянными коэффициентами k(x, z) = k:
// d^2u/dx^2 + d^2u/dz^2 + f(x, z)/k = 0
// На границах прямоугольной области (0 < x < a, 0 < z < b)

constexpr auto a = 10.0f; // см
constexpr auto b = 10.0f; // см

// Краевые условия:
// x = 0, u(0, z) = u_0,
// x = a, u(a, z) = u_0,
// z = 0, u(x, 0) = u_0,
// z = b, u(x, b) = u_0.

constexpr auto u_0 = 300.0f; // К

// В качестве функции источников можно предложить распределение вида
// f(x, z) = f_0 e^(beta (x - a/2)^2 (z - b/2)^2),
// параметры f_0, beta варьируются исходя из условия, чтобы максимум функции не превышал 3000 К.

constexpr auto f_0 = 100.0f;
constexpr auto beta = -0.001f;

constexpr auto k_0 = 2.36f; // Вт / (см * град) — теплопроводность аллюминия при комнатной температуре.

constexpr auto eps = 0.001f;

template <typename T>
constexpr T sqr(T value) {
	return value * value;
}

constexpr float f(float x, float z) {
	return f_0 * std::exp(beta * sqr(x - a / 2.0f) * sqr(z - b / 2.0f));
}

template <typename T>
using Matrix = std::vector<std::vector<T>>;

template <typename T>
void PrintMatrix(const Matrix<T>& matrix) {
	for (auto&& row : matrix) {
		for (auto&& item : row) {
			std::cout << std::setw(12) << item;
		}
		std::cout << '\n';
	}
}

void iterate() {
	constexpr auto h = 1.0f;
	constexpr auto N = static_cast<int>(a / h) + 1;
	static_assert(a == b, "If a != b use h_x, h_z, N_x and N_z");
	static_assert(N == 11);

	auto y_k = Matrix<float>(N, std::vector(N, u_0));
	auto y_k_plus_1 = y_k;
	auto y_k_max = u_0;
	auto y_k_plus_1_max = 0.0f;

	for (size_t k = 0; /* no condition */; ++k) {
		std::cout << "y^" << k << " =\n";
		PrintMatrix(y_k);

		for (size_t i = 1; i < N - 1; ++i) {
			const auto xi = static_cast<float>(i) * h;
			for (size_t j = 1; j < N - 1; ++j) {
				const auto zj = static_cast<float>(j) * h;
				y_k_plus_1[i][j] = (y_k[i + 1][j] + y_k[i - 1][j] + y_k[i][j + 1] + y_k[i][j - 1] + f(xi, zj) / k_0 * sqr(h)) / 4;
				y_k_plus_1_max = std::max(std::abs(y_k_plus_1[i][j]), y_k_plus_1_max);
			}
		}

		std::cout << "max(y^" << k << ") = " << y_k_max << std::endl;
		std::cout << "max(y^" << k + 1 << ") = " << y_k_plus_1_max << std::endl;

		if (std::abs(y_k_plus_1_max - y_k_max) <= eps) {
			break;
		}

		y_k = y_k_plus_1;

		y_k_max = y_k_plus_1_max;
		y_k_plus_1_max = 0.0f;
	}
}


int main() {
	iterate();
}
