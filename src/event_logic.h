#pragma once

#include "event.h"
#include "table_record.h"

#include <cmath>
#include <optional>

inline std::optional<EventDirection> event_direction_for_record(
	const TableRecord& rec,
	const TableRecord& obs
) {
	if (rec.s == obs.s) {
		return std::nullopt;
	}

	if (rec.s < obs.s) {
		return EventDirection::Enter;
	}

	return EventDirection::Leave;
}

inline std::optional<Real> critical_t_for_record(
	const TableRecord& rec,
	const TableRecord& obs
) {
	if (rec.s == obs.s) {
		return std::nullopt;
	}

	const Real delta_log_weight = rec.log_weight - obs.log_weight;
	const Real delta_s = static_cast<Real>(rec.s - obs.s);
	return -delta_log_weight / delta_s;
}

inline std::optional<Real> critical_p_for_record(
	const TableRecord& rec,
	const TableRecord& obs
) {
	const auto t_opt = critical_t_for_record(rec, obs);
	if (!t_opt.has_value()) {
		return std::nullopt;
	}

	const Real t = *t_opt;
	const Real exp_t = std::exp(t);
	return exp_t / (1.0L + exp_t);
}