# Command line interface

## Usage

```text
fish_sup r c n1 ... nc x1 ... xc info
```

## Meaning of the arguments

- `r`: number of rows; currently only `2` is supported
- `c`: number of columns
- `n1 ... nc`: column sums
- `x1 ... xc`: observed top-row counts
- `info`: diagnostic flag; accepted values are `true`, `false`, and `i`

The bottom row is computed internally as `n_j - x_j`.

## Example

```text
./build/fish_sup 2 3 30 40 30 20 20 10 false
```

## Output

The program prints the unconditional exact p-value together with information on where the supremum is attained.
