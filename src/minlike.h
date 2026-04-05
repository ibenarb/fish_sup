#pragma once

#include "real.h"
#include "table_record.h"

#include <cmath>
#include <stdexcept>

inline Real logit_real(Real p) {
	if (!(p > 0.0L && p < 1.0L)) {
		throw std::invalid_argument("p must satisfy 0 < p < 1");
	}

	return std::log(p / (1.0L - p));
}

inline Real minlike_log_difference(const TableRecord& rec, const TableRecord& obs, Real p) {
	return (rec.log_weight - obs.log_weight)
		 + static_cast<Real>(rec.s - obs.s) * logit_real(p);
}

inline bool is_in_minlike_set(const TableRecord& rec, const TableRecord& obs, Real p) {
	return minlike_log_difference(rec, obs, p) <= 0.0L;
}