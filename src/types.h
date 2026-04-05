#pragma once

#include <stdexcept>
#include <vector>

struct TableInput {
	int rows;
	int cols;
	std::vector<int> cells;
};

inline void validate_table_input(const TableInput& input) {
	if (input.rows < 1 || input.cols < 1) {
		throw std::invalid_argument("rows and cols must be positive");
	}

	if (static_cast<int>(input.cells.size()) != input.rows * input.cols) {
		throw std::invalid_argument("cells.size() does not match rows * cols");
	}

	for (int x : input.cells) {
		if (x < 0) {
			throw std::invalid_argument("cell entries must be nonnegative");
		}
	}
}

inline std::vector<int> row_sums(const TableInput& input) {
	validate_table_input(input);

	std::vector<int> sums(input.rows, 0);

	for (int i = 0; i < input.rows; ++i) {
		for (int j = 0; j < input.cols; ++j) {
			sums[i] += input.cells[i * input.cols + j];
		}
	}

	return sums;
}

inline std::vector<int> col_sums(const TableInput& input) {
	validate_table_input(input);

	std::vector<int> sums(input.cols, 0);

	for (int i = 0; i < input.rows; ++i) {
		for (int j = 0; j < input.cols; ++j) {
			sums[j] += input.cells[i * input.cols + j];
		}
	}

	return sums;
}

inline int total_sum(const TableInput& input) {
	validate_table_input(input);

	int total = 0;
	for (int x : input.cells) {
		total += x;
	}

	return total;
}