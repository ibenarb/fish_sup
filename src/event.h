#pragma once

#include "real.h"

#include <cstddef>

enum class EventDirection {
	Enter,
	Leave
};

struct Event {
	Real t_crit;
	Real p_crit;
	std::size_t record_index;
	EventDirection direction;
};