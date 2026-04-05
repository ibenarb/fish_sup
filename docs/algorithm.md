# Algorithm overview

fish_sup computes an exact unconditional two-sided p-value for 2 x c contingency tables under a product-binomial null model.

## Model

- Independent random variables X_j ~ Bin(n_j, p)
- Common nuisance parameter p in (0,1)
- Observed table defines the minlike acceptance region

## Core ideas

- Exact combinatorial weights are represented with GMP integers
- Tables are ordered by a minlike criterion
- The nuisance parameter is reparameterized via t = log(p / (1 - p))
- Critical values t_crit arise when the ordering of two tables changes
- These critical values are processed by an event sweep
- Between critical points, the active minlike set is constant
- The unconditional p-value is the supremum over p

## Numerical design

- Logarithms of weights are formed from sums of log binomial coefficients
- Probability evaluation is performed in a numerically stable way
- The implementation avoids lgamma-based weight construction
- Cross-platform reproducibility has been an explicit design goal
