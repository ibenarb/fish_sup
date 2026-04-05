#pragma once

#include "real.h"
#include "two_by_c.h"

#include <cmath>
#include <gmp.h>
#include <gmpxx.h>
#include <stdexcept>
#include <vector>
#pragma once

inline mpz_class binomial_mpz(unsigned long n, unsigned long k) {
    mpz_class result;
    mpz_bin_uiui(result.get_mpz_t(), n, k);
    return result;
}

inline Real log_mpz_real(const mpz_class& x) {
    if (x <= 0) {
        throw std::invalid_argument("log_mpz_real requires a positive integer");
    }

    signed long exp2 = 0;
    const double mantissa = mpz_get_d_2exp(&exp2, x.get_mpz_t());

    return ::logl(static_cast<Real>(mantissa))
         + static_cast<Real>(exp2) * ::logl(static_cast<Real>(2.0L));
}

inline mpz_class table_weight_mpz(const std::vector<int>& col_totals, const std::vector<int>& top_row) {
    if (col_totals.size() != top_row.size()) {
        throw std::invalid_argument("col_totals and top_row must have the same length");
    }

    mpz_class weight = 1;

    for (std::size_t j = 0; j < col_totals.size(); ++j) {
        const int n_j = col_totals[j];
        const int a_j = top_row[j];

        if (n_j < 0 || a_j < 0 || a_j > n_j) {
            throw std::invalid_argument("invalid column total or top-row entry");
        }

        weight *= binomial_mpz(static_cast<unsigned long>(n_j), static_cast<unsigned long>(a_j));
    }

    return weight;
}

inline mpz_class table_weight_mpz(const TwoByCData& data, const std::vector<int>& top_row) {
    return table_weight_mpz(data.col_totals, top_row);
}

inline Real log_binomial_real(unsigned long n, unsigned long k) {
    return log_mpz_real(binomial_mpz(n, k));
}

inline Real table_log_weight_real(const std::vector<int>& col_totals, const std::vector<int>& top_row) {
    if (col_totals.size() != top_row.size()) {
        throw std::invalid_argument("col_totals and top_row must have the same length");
    }

    Real value = 0.0L;

    for (std::size_t j = 0; j < col_totals.size(); ++j) {
        const int n_j = col_totals[j];
        const int a_j = top_row[j];

        if (n_j < 0 || a_j < 0 || a_j > n_j) {
            throw std::invalid_argument("invalid column total or top-row entry");
        }

        value += log_binomial_real(static_cast<unsigned long>(n_j), static_cast<unsigned long>(a_j));
    }

    return value;
}

inline Real table_log_weight_real(const TwoByCData& data, const std::vector<int>& top_row) {
    return table_log_weight_real(data.col_totals, top_row);
}