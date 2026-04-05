# fish_sup

**fish_sup** is a C++ program for computing an exact unconditional two-sided test for 2 x c contingency tables under a product-binomial model.

## Overview

The program implements an exact unconditional two-sided test of Barnard type for 2 x c tables under independent binomial sampling with a common nuisance parameter.

Core characteristics:

- exact
- unconditional
- two-sided
- minlike-based ordering
- supremum over p in (0,1)
- event-sweep evaluation
- GMP-backed exact combinatorial weights

## Quick start

Build:

```text
cmake -S . -B build && cmake --build build
```

Run the reference example:

```text
./build/fish_sup 2 3 30 40 30 20 20 10 false
```

## Reference result

For the command above, the current implementation returns:

```text
unconditional test result:
  p_value = 0.0863981038906587926911
  attained_as = left_limit
  at p = 0.501295727632370348025
  at t = 0.00518292213171590368517
  block index = 9727
```

## Documentation

- [Algorithm overview](docs/algorithm.md)
- [Numerical notes](docs/numerics.md)
- [Reference example](examples/reference_2x3.md)
- [CLI documentation](docs/cli.md)

## Build requirements

- CMake >= 3.28
- C++20 compiler
- GMP
- GMP C++ bindings (`gmpxx`)

## Status

The repository currently contains a working command-line implementation together with initial documentation and a reproducible reference example.

## License

This project is released under the MIT License. See [LICENSE](LICENSE).
