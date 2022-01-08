#include <chrono>
#include <iostream>

#include <io/print_matrix.hpp>
#include <model/finite_difference_approximation.hpp>
#include <model/probability.hpp>


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

constexpr auto eps = 0.0001f;


void iterate() {
	constexpr auto h = 1.0f;
	constexpr auto N = static_cast<int>(a / h) + 1;
	static_assert(a == b, "If a != b use h_x, h_z, N_x and N_z");
	static_assert(N == 11);

	constexpr Parameters parameters = {
		.a = a,
		.b = b,

		.N_x = N,
		.N_z = N,

		.tau = 1.0f,
		.eps = eps,

		.u_0 = u_0,

		.f_0 = f_0,
		.beta = beta,

		.k = k_0,

		.N_p = 1000,
	};

	auto t1 = std::chrono::high_resolution_clock::now();
	auto dependency = FDA::solve(parameters);
	auto t2 = std::chrono::high_resolution_clock::now();

	PrintMatrix(dependency.y.back());
	std::cout << dependency.t.back() << std::endl;
	std::cout << std::chrono::duration<double, std::milli>{t2 - t1}.count() << std::endl;

	auto diff = dependency.y.back();

	auto t3 = std::chrono::high_resolution_clock::now();
	dependency = P::solve(parameters);
	auto t4 = std::chrono::high_resolution_clock::now();

	PrintMatrix(dependency.y.back());
	std::cout << dependency.t.back() << std::endl;
	std::cout << std::chrono::duration<double, std::milli>{t4 - t3}.count() << std::endl;

	for (Size i = 0; i < diff.size(); ++i) {
		for (Size j = 0; j < diff[i].size(); ++j) {
			diff[i][j] -= dependency.y.back()[i][j];
		}
	}
	PrintMatrix(diff);
}


int main() {
	iterate();
}
