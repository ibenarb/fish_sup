#pragma once

#include "two_by_c.h"

#include <functional>
#include <vector>

inline void enumerate_top_rows_recursive(
	const std::vector<int>& col_totals,
	int j,
	std::vector<int>& current,
	const std::function<void(const std::vector<int>&)>& callback
) {
	if (j == static_cast<int>(col_totals.size())) {
		callback(current);
		return;
	}

	for (int a = 0; a <= col_totals[j]; ++a) {
		current[j] = a;
		enumerate_top_rows_recursive(col_totals, j + 1, current, callback);
	}
}

inline void enumerate_top_rows(
	const TwoByCData& data,
	const std::function<void(const std::vector<int>&)>& callback
) {
	std::vector<int> current(data.cols, 0);
	enumerate_top_rows_recursive(data.col_totals, 0, current, callback);
}