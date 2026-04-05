#pragma once

#include "real.h"

#include <cmath>

inline Real logsumexp_pair(Real a, Real b) {
	if (a == -INFINITY) {
		return b;
	}
	if (b == -INFINITY) {
		return a;
	}
	const Real m = (a > b) ? a : b;
	return m + std::log(std::exp(a - m) + std::exp(b - m));
}