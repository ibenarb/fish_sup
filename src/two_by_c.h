#pragma once

#include "types.h"

#include <numeric>
#include <stdexcept>
#include <vector>

struct TwoByCData {
	int cols;
	std::vector<int> top;
	std::vector<int> bottom;
	std::vector<int> col_totals;
	int top_sum;
	int bottom_sum;
	int total;
};

inline TwoByCData make_two_by_c_data(const TableInput& input) {
	validate_table_input(input);

	if (input.rows != 2) {
		throw std::invalid_argument("table must have exactly 2 rows");
	}

	if (input.cols < 1) {
		throw std::invalid_argument("table must have at least 1 column");
	}

	TwoByCData data{};
	data.cols = input.cols;
	data.top.resize(input.cols);
	data.bottom.resize(input.cols);
	data.col_totals = col_sums(input);

	for (int j = 0; j < input.cols; ++j) {
		data.top[j] = input.cells[j];
		data.bottom[j] = input.cells[input.cols + j];
	}

	data.top_sum = 0;
	data.bottom_sum = 0;
	for (int x : data.top) {
		data.top_sum += x;
	}
	for (int x : data.bottom) {
		data.bottom_sum += x;
	}

	data.total = data.top_sum + data.bottom_sum;
	return data;
}

inline const std::vector<int>& observed_top_row(const TwoByCData& data) {
	return data.top;
}

inline int observed_success_sum(const TwoByCData& data) {
	return data.top_sum;
}