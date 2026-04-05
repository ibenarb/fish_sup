#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "event.h"
#include "event_logic.h"
#include "table_record.h"

struct EventBlock {
    long double t_crit = 0.0L;
    std::vector<std::size_t> event_indices;
    std::size_t enter_count = 0;
    std::size_t leave_count = 0;
};

struct EventCollection {
    std::vector<Event> events;
    std::vector<EventBlock> blocks;
};

struct TGroupingConfig {
    long double abs_tol = 1e-18L;
    long double rel_tol = 1e-15L;
};

inline bool nearly_equal_t(long double a, long double b, const TGroupingConfig& cfg) {
    const long double diff = std::fabs(a - b);
    const long double scale = std::max(std::fabs(a), std::fabs(b));
    return diff <= cfg.abs_tol + cfg.rel_tol * scale;
}

inline std::vector<Event> build_all_events(
    const std::vector<TableRecord>& records,
    const TableRecord& obs
) {
    std::vector<Event> events;
    events.reserve(records.size());

    for (std::size_t i = 0; i < records.size(); ++i) {
        const auto dir_opt = event_direction_for_record(records[i], obs);
        if (!dir_opt.has_value()) {
            continue;
        }

        const auto t_opt = critical_t_for_record(records[i], obs);
        const auto p_opt = critical_p_for_record(records[i], obs);

        if (!t_opt.has_value() || !p_opt.has_value()) {
            continue;
        }

        events.push_back(Event{
            *t_opt,
            *p_opt,
            i,
            *dir_opt
        });
    }

    std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
        if (a.t_crit < b.t_crit) {
            return true;
        }
        if (b.t_crit < a.t_crit) {
            return false;
        }
        if (a.direction != b.direction) {
            return a.direction == EventDirection::Enter;
        }
        return a.record_index < b.record_index;
    });

    return events;
}

inline std::vector<EventBlock> build_event_blocks(
    const std::vector<Event>& events,
    const std::vector<TableRecord>& records,
    const TGroupingConfig& cfg
) {
    (void)records;

    std::vector<EventBlock> blocks;

    if (events.empty()) {
        return blocks;
    }

    EventBlock current;
    current.t_crit = events[0].t_crit;
    current.event_indices.push_back(0);
    if (events[0].direction == EventDirection::Enter) {
        current.enter_count = 1;
    } else {
        current.leave_count = 1;
    }

    for (std::size_t i = 1; i < events.size(); ++i) {
        const long double prev_t = events[i - 1].t_crit;
        const long double curr_t = events[i].t_crit;

        if (nearly_equal_t(prev_t, curr_t, cfg)) {
            current.event_indices.push_back(i);
            if (events[i].direction == EventDirection::Enter) {
                ++current.enter_count;
            } else {
                ++current.leave_count;
            }
        } else {
            blocks.push_back(current);

            current = EventBlock{};
            current.t_crit = events[i].t_crit;
            current.event_indices.push_back(i);
            if (events[i].direction == EventDirection::Enter) {
                current.enter_count = 1;
            } else {
                current.leave_count = 1;
            }
        }
    }

    blocks.push_back(current);
    return blocks;
}

inline EventCollection build_event_collection(
    const std::vector<TableRecord>& records,
    const TableRecord& obs,
    const TGroupingConfig& cfg
) {
    EventCollection result;
    result.events = build_all_events(records, obs);
    result.blocks = build_event_blocks(result.events, records, cfg);
    return result;
}