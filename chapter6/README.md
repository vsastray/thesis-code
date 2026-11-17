# Chapter 6 — Figure 6.2

Capacity calculations for two identical amplitude-damping channels, comparing `Free` and `FreeDef`.

Requires MATLAB, CVX, QETLAB, and quantinf; see [setup](../README.md#requirements).

## Run

For one parameter point, set this folder as MATLAB's current folder and enter:

```matlab
clear;
damp_rate = 0.3;
epsilon = 0.02;
classical_capacity_ico;
```

Inspect `capacity_ico`, `capacity_def`, `cvx_status_ico`, and `cvx_status_def`. The objective values are the parameter `m`; the figure uses `log2(floor(m))` bits.

Run [collect_capacity_data.m](collect_capacity_data.m) for damping rates `0:0.01:0.5` at error tolerances `0` and `0.02`. It writes `capacity_data.mat` and `capacity_with_err_data.mat`, overwriting existing files. The files contain objective values only; solver statuses are not saved.

## Numerical notes

Three parameter points were checked; the full sweep has not been rerun. The two zero-error checks returned `Inaccurate/Solved`; the bounded-error check returned `Solved`.

Integer thresholds require care. The original plot set the zero-error `Free` endpoint at zero damping to 4 bits because the saved objective was slightly below 16.

## Related paper

Xuanqiang Zhao, Benchi Zhao, Cyril Branciard, and Giulio Chiribella. [*The communication power of indefinite causal order*](https://arxiv.org/abs/2510.08507) (2025).
