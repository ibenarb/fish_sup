#pragma once

#include "enumeration.h"
#include "table_record.h"

#include <stdexcept>
#include <vector>

inline std::size_t number_of_top_rows(const TwoByCData& data) {
	std::size_t count = 1;

	for (int n_j : data.col_totals) {
		if (n_j < 0) {
			throw std::invalid_argument("column totals must be nonnegative");
		}
		count *= static_cast<std::size_t>(n_j + 1);
	}

	return count;
}

inline std::vector<TableRecord> build_all_table_records(const TwoByCData& data) {
	std::vector<TableRecord> records;
	records.reserve(number_of_top_rows(data));

	enumerate_top_rows(data, [&](const std::vector<int>& top_row) {
		records.push_back(make_table_record(data, top_row));
	});

	return records;
}