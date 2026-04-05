#pragma once

#include "combinatorics.h"
#include "event_collection.h"
#include "minlike.h"
#include "real.h"
#include "table_record.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

enum class ThresholdInclusionMode {
    OffThreshold,
    OnThreshold
};

struct SweepWeightState {
    std::vector<mpz_class> weights_by_s;
};

inline Real record_log_probability(
    const TableRecord& rec,
    int total_n,
    Real p
) {
    if (!(p > 0.0L && p < 1.0L)) {
        throw std::invalid_argument("p must satisfy 0 < p < 1");
    }

    const Real log_p = ::logl(p);
    const Real log_q = ::log1pl(-p);

    return rec.log_weight
         + static_cast<Real>(rec.s) * log_p
         + static_cast<Real>(total_n - rec.s) * log_q;
}

inline SweepWeightState initial_minlike_weight_state_left_of_all_thresholds(
    const std::vector<TableRecord>& records,
    const TableRecord& obs,
    int total_n
) {
    SweepWeightState state;
    state.weights_by_s.assign(static_cast<std::size_t>(total_n + 1), 0);

    for (const TableRecord& rec : records) {
        if (rec.s > obs.s || (rec.s == obs.s && rec.log_weight <= obs.log_weight)) {
            state.weights_by_s[static_cast<std::size_t>(rec.s)] += rec.weight;
        }
    }

    return state;
}

inline void apply_event_block_strictly_left_to_weight_state(
    SweepWeightState& state,
    const EventBlock& block,
    const std::vector<Event>& events,
    const std::vector<TableRecord>& records
) {
    for (std::size_t event_index : block.event_indices) {
        const Event& ev = events[event_index];
        const TableRecord& rec = records[ev.record_index];
        mpz_class& slot = state.weights_by_s[static_cast<std::size_t>(rec.s)];

        if (ev.direction == EventDirection::Enter) {
            slot += rec.weight;
        } else {
            if (slot < rec.weight) {
                throw std::runtime_error("negative weight state in leave update");
            }
            slot -= rec.weight;
        }
    }
}

inline void apply_event_block_on_threshold_to_weight_state(
    SweepWeightState& state,
    const EventBlock& block,
    const std::vector<Event>& events,
    const std::vector<TableRecord>& records
) {
    for (std::size_t event_index : block.event_indices) {
        const Event& ev = events[event_index];
        if (ev.direction != EventDirection::Enter) {
            continue;
        }

        const TableRecord& rec = records[ev.record_index];
        state.weights_by_s[static_cast<std::size_t>(rec.s)] += rec.weight;
    }
}

inline Real weight_array_probability(
    const std::vector<mpz_class>& weights,
    int total_n,
    Real p
) {
    if (!(p > 0.0L && p < 1.0L)) {
        throw std::invalid_argument("p must satisfy 0 < p < 1");
    }

    if (static_cast<int>(weights.size()) != total_n + 1) {
        throw std::invalid_argument("weights.size() must be total_n + 1");
    }

    const Real log_p = ::logl(p);
    const Real log_q = ::log1pl(-p);

    bool has_term = false;
    Real max_log = 0.0L;

    for (int s = 0; s <= total_n; ++s) {
        const mpz_class& w = weights[static_cast<std::size_t>(s)];
        if (w == 0) {
            continue;
        }

        const Real log_term =
            log_mpz_real(w)
            + static_cast<Real>(s) * log_p
            + static_cast<Real>(total_n - s) * log_q;

        if (!has_term || log_term > max_log) {
            max_log = log_term;
            has_term = true;
        }
    }

    if (!has_term) {
        return 0.0L;
    }

    Real scaled_sum = 0.0L;

    for (int s = 0; s <= total_n; ++s) {
        const mpz_class& w = weights[static_cast<std::size_t>(s)];
        if (w == 0) {
            continue;
        }

        const Real log_term =
            log_mpz_real(w)
            + static_cast<Real>(s) * log_p
            + static_cast<Real>(total_n - s) * log_q;

        scaled_sum += ::expl(log_term - max_log);
    }

    return ::expl(max_log + ::logl(scaled_sum));
}

inline Real minlike_p_value_event_sweep_exact(
    const std::vector<TableRecord>& records,
    const TableRecord& obs,
    const EventCollection& ec,
    int total_n,
    Real p,
    const TGroupingConfig& cfg,
    ThresholdInclusionMode mode
) {
    if (!(p > 0.0L && p < 1.0L)) {
        throw std::invalid_argument("p must satisfy 0 < p < 1");
    }

    SweepWeightState state =
        initial_minlike_weight_state_left_of_all_thresholds(records, obs, total_n);

    const Real t = logit_real(p);

    for (const EventBlock& block : ec.blocks) {
        const bool is_equal = nearly_equal_t(block.t_crit, t, cfg);
        const bool is_strict_left = (block.t_crit < t) && !is_equal;

        if (is_strict_left) {
            apply_event_block_strictly_left_to_weight_state(state, block, ec.events, records);
            continue;
        }

        if (is_equal && mode == ThresholdInclusionMode::OnThreshold) {
            apply_event_block_on_threshold_to_weight_state(state, block, ec.events, records);
        }

        break;
    }

    return weight_array_probability(state.weights_by_s, total_n, p);
}

inline Real minlike_p_value_event_sweep_exact_auto(
    const std::vector<TableRecord>& records,
    const TableRecord& obs,
    const EventCollection& ec,
    int total_n,
    Real p,
    const TGroupingConfig& cfg
) {
    if (!(p > 0.0L && p < 1.0L)) {
        throw std::invalid_argument("p must satisfy 0 < p < 1");
    }

    const Real t = logit_real(p);
    ThresholdInclusionMode mode = ThresholdInclusionMode::OffThreshold;

    for (const EventBlock& block : ec.blocks) {
        if (nearly_equal_t(block.t_crit, t, cfg)) {
            mode = ThresholdInclusionMode::OnThreshold;
            break;
        }
        if (t < block.t_crit && !nearly_equal_t(block.t_crit, t, cfg)) {
            break;
        }
    }

    return minlike_p_value_event_sweep_exact(records, obs, ec, total_n, p, cfg, mode);
}

inline std::vector<mpz_class> minlike_weight_array_bruteforce(
    const std::vector<TableRecord>& records,
    const TableRecord& obs,
    int total_n,
    Real p
) {
    std::vector<mpz_class> weights(static_cast<std::size_t>(total_n + 1), 0);

    for (const TableRecord& rec : records) {
        if (is_in_minlike_set(rec, obs, p)) {
            weights[static_cast<std::size_t>(rec.s)] += rec.weight;
        }
    }

    return weights;
}

inline Real minlike_p_value_bruteforce(
    const std::vector<TableRecord>& records,
    const TableRecord& obs,
    int total_n,
    Real p
) {
    const std::vector<mpz_class> weights =
        minlike_weight_array_bruteforce(records, obs, total_n, p);

    return weight_array_probability(weights, total_n, p);
}