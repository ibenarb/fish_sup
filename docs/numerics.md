# Numerical notes

The implementation of fish_sup was designed with numerical reproducibility in mind.

## Main points

- Exact combinatorial weights are stored as GMP integers (mpz_class)
- Log weights are evaluated as sums of log binomial coefficients
- No lgamma-based reconstruction of binomial weights is used
- Probability aggregation is carried out in a numerically stable way
- The code has been checked for consistent results across Linux and Windows builds

## Reproducibility

For the documented 2 x 3 reference example, the current implementation returns:

```text
p_value = 0.0863981038906587926911
attained_as = left_limit
```

Small differences in final printed digits may occur if compiler options, floating-point types, or library behavior are changed. The current build has been tuned to avoid unnecessary platform-dependent drift.
