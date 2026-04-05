# fish_sup

**fish_sup** is a C++ program for computing an exact unconditional two-sided test for 2 x c contingency tables under a product-binomial model.

## Quick start

Build:

cmake -S . -B build && cmake --build build

Run example:

./build/fish_sup 2 3 30 40 30 20 20 10 false

## Reference result

p_value = 0.0863981038906587926911

## Status

Exact, unconditional, minlike-based test with event sweep and GMP-backed arithmetic.
