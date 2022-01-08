#include <io/print_matrix.hpp>

#include <iostream>
#include <iomanip>


void PrintMatrix(const Matrix& matrix, int item_width) {
	for (auto&& row : matrix) {
		for (auto&& item : row) {
			std::cout << std::setw(item_width) << item;
		}
		std::cout << '\n';
	}
}

void PrintMatrixW(const Matrix& matrix, int item_width) {
	for (auto&& row : matrix) {
		for (auto&& item : row) {
			std::cout << std::setw(item_width) << std::fixed << std::setprecision(4) << item;
		}
		std::cout << '\n';
	}
}
