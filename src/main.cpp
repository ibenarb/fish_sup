#include "event_collection.h"
#include "minlike.h"
#include "probability.h"
#include "table_collection.h"
#include "two_by_c.h"
#include "types.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

struct ExactSweepProfilePoint {
    Real t_crit = 0.0L;
    Real p_crit = 0.0L;
    Real value_at_threshold = 0.0L;
    std::vector<mpz_class> weights_strictly_right;
};

struct ExactSweepProfile {
    std::vector<mpz_class> weights_left_of_all_thresholds;
    std::vector<ExactSweepProfilePoint> points;
};

struct UnconditionalSupremumResult {
    std::string attained_as;
    Real value = -1.0L;
    std::size_t index = 0;
    Real t = 0.0L;
    Real p = 0.0L;
};

struct CliInput {
    int rows = 0;
    int cols = 0;
    std::vector<int> column_sums;
    std::vector<int> top_row;
    bool diagnostics = false;
};

inline const char* build_version() {
    return __DATE__ " " __TIME__;
}

inline Real logistic_from_t(Real t) {
    if (t >= static_cast<Real>(0)) {
        const Real e = ::expl(-t);
        return static_cast<Real>(1) / (static_cast<Real>(1) + e);
    } else {
        const Real e = ::expl(t);
        return e / (static_cast<Real>(1) + e);
    }
}

inline ExactSweepProfile build_exact_sweep_profile(
    const std::vector<TableRecord>& records,
    const TableRecord& obs,
    const EventCollection& ec,
    int total_n
) {
    ExactSweepProfile profile;

    profile.weights_left_of_all_thresholds =
        initial_minlike_weight_state_left_of_all_thresholds(records, obs, total_n).weights_by_s;

    SweepWeightState state_left;
    state_left.weights_by_s = profile.weights_left_of_all_thresholds;

    profile.points.reserve(ec.blocks.size());

    for (const EventBlock& block : ec.blocks) {
        ExactSweepProfilePoint point;
        point.t_crit = block.t_crit;
        point.p_crit = logistic_from_t(block.t_crit);

        SweepWeightState state_at_threshold = state_left;
        apply_event_block_on_threshold_to_weight_state(state_at_threshold, block, ec.events, records);
        point.value_at_threshold =
            weight_array_probability(state_at_threshold.weights_by_s, total_n, point.p_crit);

        SweepWeightState state_right = state_left;
        apply_event_block_strictly_left_to_weight_state(state_right, block, ec.events, records);
        point.weights_strictly_right = state_right.weights_by_s;

        profile.points.push_back(point);
        state_left = state_right;
    }

    return profile;
}

inline Real midpoint(Real a, Real b) {
    return 0.5L * (a + b);
}

inline UnconditionalSupremumResult find_unconditional_supremum_on_profile(
    const ExactSweepProfile& profile,
    int total_n
) {
    UnconditionalSupremumResult best_left{"left_limit", -1.0L, 0, 0.0L, 0.0L};
    UnconditionalSupremumResult best_threshold{"threshold", -1.0L, 0, 0.0L, 0.0L};
    UnconditionalSupremumResult best_right{"right_limit", -1.0L, 0, 0.0L, 0.0L};

    for (std::size_t i = 0; i < profile.points.size(); ++i) {
        const Real p = profile.points[i].p_crit;
        const Real t = profile.points[i].t_crit;

        Real left_limit = 0.0L;
        if (i == 0) {
            left_limit = weight_array_probability(profile.weights_left_of_all_thresholds, total_n, p);
        } else {
            left_limit = weight_array_probability(profile.points[i - 1].weights_strictly_right, total_n, p);
        }

        const Real threshold_value = profile.points[i].value_at_threshold;
        const Real right_limit =
            weight_array_probability(profile.points[i].weights_strictly_right, total_n, p);

        if (left_limit > best_left.value) {
            best_left.value = left_limit;
            best_left.index = i;
            best_left.t = t;
            best_left.p = p;
        }

        if (threshold_value > best_threshold.value) {
            best_threshold.value = threshold_value;
            best_threshold.index = i;
            best_threshold.t = t;
            best_threshold.p = p;
        }

        if (right_limit > best_right.value) {
            best_right.value = right_limit;
            best_right.index = i;
            best_right.t = t;
            best_right.p = p;
        }
    }

    UnconditionalSupremumResult best = best_left;
    if (best_threshold.value > best.value) {
        best = best_threshold;
    }
    if (best_right.value > best.value) {
        best = best_right;
    }

    return best;
}

inline UnconditionalSupremumResult unconditional_supremum_result(
    const ExactSweepProfile& profile,
    int total_n
) {
    return find_unconditional_supremum_on_profile(profile, total_n);
}

inline UnconditionalSupremumResult unconditional_supremum_result(
    const std::vector<TableRecord>& records,
    const TableRecord& obs,
    const EventCollection& ec,
    int total_n
) {
    const ExactSweepProfile profile = build_exact_sweep_profile(records, obs, ec, total_n);
    return unconditional_supremum_result(profile, total_n);
}

inline void print_unconditional_test_result(const UnconditionalSupremumResult& sup) {
    std::cout << "unconditional test result:\n";
    std::cout << "  p_value = " << sup.value << '\n';
    std::cout << "  attained_as = " << sup.attained_as << '\n';
    std::cout << "  at p = " << sup.p << '\n';
    std::cout << "  at t = " << sup.t << '\n';
    std::cout << "  block index = " << sup.index << '\n';
}

inline void run_reference_diagnostics(
    const std::vector<TableRecord>& records,
    const TableRecord& obs,
    const EventCollection& ec,
    int total_n,
    const TGroupingConfig& cfg
) {
    const ExactSweepProfile profile = build_exact_sweep_profile(records, obs, ec, total_n);
    const UnconditionalSupremumResult sup_from_profile = unconditional_supremum_result(profile, total_n);
    const UnconditionalSupremumResult sup_direct = unconditional_supremum_result(records, obs, ec, total_n);

    const std::vector<Real> test_p{
        1e-30L,
        1e-18L,
        0.01L,
        0.1L,
        0.5L,
        0.9L,
        0.99L
    };

    std::cout << "record count = " << records.size() << '\n';
    std::cout << "event count = " << ec.events.size() << '\n';
    std::cout << "block count = " << ec.blocks.size() << '\n';
    std::cout << "profile point count = " << profile.points.size() << '\n';

    Real max_profile_threshold_diff = 0.0L;
    std::size_t max_profile_threshold_diff_index = 0;

    for (std::size_t i = 0; i < profile.points.size(); ++i) {
        const Real p = profile.points[i].p_crit;
        const Real sweep_at_p =
            minlike_p_value_event_sweep_exact_auto(records, obs, ec, total_n, p, cfg);
        const Real diff = std::fabs(profile.points[i].value_at_threshold - sweep_at_p);

        if (diff > max_profile_threshold_diff) {
            max_profile_threshold_diff = diff;
            max_profile_threshold_diff_index = i;
        }
    }

    std::cout << "max profile threshold diff = " << max_profile_threshold_diff << '\n';
    std::cout << "max profile threshold diff index = " << max_profile_threshold_diff_index << '\n';

    Real max_profile_right_diff = 0.0L;
    std::size_t max_profile_right_diff_index = 0;
    std::size_t checked_open_intervals = 0;

    Real best_midpoint_value = -1.0L;
    std::size_t best_midpoint_index = 0;
    Real best_midpoint_t = 0.0L;
    Real best_midpoint_p = 0.0L;

    for (std::size_t i = 0; i + 1 < profile.points.size(); ++i) {
        const Real t_left = profile.points[i].t_crit;
        const Real t_right = profile.points[i + 1].t_crit;

        if (nearly_equal_t(t_left, t_right, cfg)) {
            continue;
        }

        const Real t_mid = midpoint(t_left, t_right);
        const Real p_mid = logistic_from_t(t_mid);

        const Real sweep_mid =
            minlike_p_value_event_sweep_exact_auto(records, obs, ec, total_n, p_mid, cfg);

        const Real state_mid =
            weight_array_probability(profile.points[i].weights_strictly_right, total_n, p_mid);

        const Real diff = std::fabs(sweep_mid - state_mid);

        if (diff > max_profile_right_diff) {
            max_profile_right_diff = diff;
            max_profile_right_diff_index = i;
        }

        if (sweep_mid > best_midpoint_value) {
            best_midpoint_value = sweep_mid;
            best_midpoint_index = i;
            best_midpoint_t = t_mid;
            best_midpoint_p = p_mid;
        }

        ++checked_open_intervals;
    }

    std::cout << "checked open intervals = " << checked_open_intervals << '\n';
    std::cout << "max profile right-state diff = " << max_profile_right_diff << '\n';
    std::cout << "max profile right-state diff index = " << max_profile_right_diff_index << '\n';

    std::cout << "best midpoint value = " << best_midpoint_value << '\n';
    std::cout << "best midpoint interval index = " << best_midpoint_index << '\n';
    std::cout << "best midpoint t = " << best_midpoint_t << '\n';
    std::cout << "best midpoint p = " << best_midpoint_p << '\n';

    std::cout << "sup_from_profile attained_as = " << sup_from_profile.attained_as << '\n';
    std::cout << "sup_from_profile value = " << sup_from_profile.value << '\n';
    std::cout << "sup_from_profile index = " << sup_from_profile.index << '\n';
    std::cout << "sup_from_profile t = " << sup_from_profile.t << '\n';
    std::cout << "sup_from_profile p = " << sup_from_profile.p << '\n';

    std::cout << "sup_direct attained_as = " << sup_direct.attained_as << '\n';
    std::cout << "sup_direct value = " << sup_direct.value << '\n';
    std::cout << "sup_direct index = " << sup_direct.index << '\n';
    std::cout << "sup_direct t = " << sup_direct.t << '\n';
    std::cout << "sup_direct p = " << sup_direct.p << '\n';

    const Real sup_value_diff = std::fabs(sup_from_profile.value - sup_direct.value);
    const Real sup_t_diff = std::fabs(sup_from_profile.t - sup_direct.t);
    const Real sup_p_diff = std::fabs(sup_from_profile.p - sup_direct.p);
    const bool sup_index_equal = (sup_from_profile.index == sup_direct.index);
    const bool sup_kind_equal = (sup_from_profile.attained_as == sup_direct.attained_as);

    std::cout << "sup value diff = " << sup_value_diff << '\n';
    std::cout << "sup t diff = " << sup_t_diff << '\n';
    std::cout << "sup p diff = " << sup_p_diff << '\n';
    std::cout << "sup index equal = " << sup_index_equal << '\n';
    std::cout << "sup kind equal = " << sup_kind_equal << '\n';

    std::cout << "unconditional p value = " << sup_direct.value << '\n';
    std::cout << "unconditional p diff to sup_direct = " << 0.0L << '\n';

    const Real sweep_at_sup =
        minlike_p_value_event_sweep_exact_auto(records, obs, ec, total_n, sup_direct.p, cfg);
    const Real sup_diff = std::fabs(sup_direct.value - sweep_at_sup);

    std::cout << "sweep at sup p = " << sweep_at_sup << '\n';
    std::cout << "sup abs diff = " << sup_diff << '\n';

    for (Real p : test_p) {
        const Real sweep = minlike_p_value_event_sweep_exact_auto(records, obs, ec, total_n, p, cfg);
        const Real brute = minlike_p_value_bruteforce(records, obs, total_n, p);

        std::cout << "p = " << p << '\n';
        std::cout << "  sweep_exact = " << sweep << '\n';
        std::cout << "  brute = " << brute << '\n';
        std::cout << "  abs diff = " << std::fabs(sweep - brute) << '\n';
    }
}

inline void print_help(const char* program_name) {
    std::cout << "Usage:\n";
    std::cout << "  " << program_name << " r c n1 ... nc x1 ... xc info\n";
    std::cout << '\n';
    std::cout << "Meaning:\n";
    std::cout << "  r       number of rows (currently must be 2)\n";
    std::cout << "  c       number of columns (positive integer)\n";
    std::cout << "  n1..nc  column sums\n";
    std::cout << "  x1..xc  observed top-row counts\n";
    std::cout << "  info    one of: true, false, i\n";
    std::cout << '\n';
    std::cout << "The bottom row is computed internally as nj - xj.\n";
    std::cout << '\n';
    std::cout << "Examples:\n";
    std::cout << "  " << program_name << " 2 3 30 40 30 20 20 10 false\n";
    std::cout << "  " << program_name << " 2 3 30 40 30 20 20 10 true\n";
    std::cout << '\n';
    std::cout << "Help:\n";
    std::cout << "  " << program_name << " ?\n";
    std::cout << "  " << program_name << " h\n";
}

inline int parse_int_strict(const std::string& text, const std::string& name) {
    std::size_t pos = 0;
    const long long value = std::stoll(text, &pos);
    if (pos != text.size()) {
        throw std::runtime_error("invalid integer for " + name + ": " + text);
    }
    if (value < static_cast<long long>(std::numeric_limits<int>::min()) ||
        value > static_cast<long long>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("integer out of range for " + name + ": " + text);
    }
    return static_cast<int>(value);
}

inline bool parse_info_flag(const std::string& text) {
    if (text == "i" || text == "true") {
        return true;
    }
    if (text == "false") {
        return false;
    }
    throw std::runtime_error("info must be one of: true, false, i");
}

inline CliInput parse_cli_input(int argc, char* argv[]) {
    if (argc == 2) {
        const std::string arg = argv[1];
        if (arg == "?" || arg == "h") {
            print_help(argv[0]);
            std::exit(0);
        }
    }

    if (argc < 5) {
        throw std::runtime_error("too few arguments");
    }

    CliInput input;
    input.rows = parse_int_strict(argv[1], "r");
    input.cols = parse_int_strict(argv[2], "c");

    if (input.rows != 2) {
        throw std::runtime_error("currently only r = 2 is supported");
    }
    if (input.cols <= 0) {
        throw std::runtime_error("c must be positive");
    }

    const int expected_argc = 2 * input.cols + 4;
    if (argc != expected_argc) {
        throw std::runtime_error(
            "wrong number of arguments: expected " + std::to_string(expected_argc - 1) +
            " user arguments after program name"
        );
    }

    input.column_sums.reserve(input.cols);
    input.top_row.reserve(input.cols);

    for (int j = 0; j < input.cols; ++j) {
        const int n_j = parse_int_strict(argv[3 + j], "n" + std::to_string(j + 1));
        if (n_j < 0) {
            throw std::runtime_error("column sums must be nonnegative");
        }
        input.column_sums.push_back(n_j);
    }

    for (int j = 0; j < input.cols; ++j) {
        const int x_j = parse_int_strict(argv[3 + input.cols + j], "x" + std::to_string(j + 1));
        if (x_j < 0) {
            throw std::runtime_error("top-row counts must be nonnegative");
        }
        if (x_j > input.column_sums[j]) {
            throw std::runtime_error("each xj must satisfy 0 <= xj <= nj");
        }
        input.top_row.push_back(x_j);
    }

    input.diagnostics = parse_info_flag(argv[3 + 2 * input.cols]);
    return input;
}

inline TableInput make_table_input_from_cli(const CliInput& cli) {
    std::vector<int> cells;
    cells.reserve(2 * cli.cols);

    for (int j = 0; j < cli.cols; ++j) {
        cells.push_back(cli.top_row[j]);
    }
    for (int j = 0; j < cli.cols; ++j) {
        cells.push_back(cli.column_sums[j] - cli.top_row[j]);
    }

    return TableInput{2, cli.cols, cells};
}

int main(int argc, char* argv[]) {
    try {
        std::cout << "build_version = " << build_version() << '\n';
        std::cout << "sizeof(Real) = " << sizeof(Real) << '\n';
        std::cout << "digits10(Real) = " << std::numeric_limits<Real>::digits10 << '\n';
        std::cout << "max_digits10(Real) = " << std::numeric_limits<Real>::max_digits10 << '\n';
        std::cout << "is_iec559(Real) = " << std::numeric_limits<Real>::is_iec559 << '\n';

        const CliInput cli = parse_cli_input(argc, argv);
        const TableInput input = make_table_input_from_cli(cli);

        TwoByCData data = make_two_by_c_data(input);
        std::vector<TableRecord> records = build_all_table_records(data);
        TableRecord obs = make_table_record(data, observed_top_row(data));

        const int total_n = data.total;

        TGroupingConfig cfg{};
        EventCollection ec = build_event_collection(records, obs, cfg);

        const UnconditionalSupremumResult sup =
            unconditional_supremum_result(records, obs, ec, total_n);

        std::cout << std::setprecision(21);
        print_unconditional_test_result(sup);

        if (cli.diagnostics) {
            run_reference_diagnostics(records, obs, ec, total_n, cfg);
        }

        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << '\n';
        print_help(argv[0]);
        return 1;
    }
}