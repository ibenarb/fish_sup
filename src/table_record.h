#pragma once

#include "combinatorics.h"
#include "real.h"

#include <vector>

struct TableRecord {
	std::vector<int> top_row;
	int s;
	mpz_class weight;
	Real log_weight;
};

inline TableRecord make_table_record(const std::vector<int>& col_totals, const std::vector<int>& top_row) {
	TableRecord rec{};
	rec.top_row = top_row;
	rec.s = 0;

	for (int x : top_row) {
		rec.s += x;
	}

	rec.weight = table_weight_mpz(col_totals, top_row);
	rec.log_weight = table_log_weight_real(col_totals, top_row);

	return rec;
}

inline TableRecord make_table_record(const TwoByCData& data, const std::vector<int>& top_row) {
	return make_table_record(data.col_totals, top_row);
}