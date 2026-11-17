# Chapter 4 — Figures 4.2–4.4

Requires MATLAB, CVX, QETLAB, and quantinf; see [setup](../README.md#requirements).

## Run

Set this folder as MATLAB's current folder. In each script, keep one `noise_channel` assignment active, set the parameters below, and run the file. Keep the other dimension settings at their defaults.

| Figure | Script | Parameter values |
| --- | --- | --- |
| 4.2 | [shadow_comm_sample_complexity.m](shadow_comm_sample_complexity.m) | `dim_target`: `2`, `3`, `4` |
| 4.3 | [noise_simulation_sample_complexity.m](noise_simulation_sample_complexity.m) | `dim_id`: `2`, `3` |
| 4.4, left | [shadow_comm_min_err.m](shadow_comm_min_err.m) | `cost_budget`: `0.90`, `0.95`, …, `1.20` |
| 4.4, right | [noise_simulation_min_err.m](noise_simulation_min_err.m) | `cost_budget`: `0.90`, `0.95`, …, `1.20` |

Run each parameter value separately for depolarizing, dephasing, and amplitude-damping noise. The thesis parameter is `p = 0.9`; amplitude damping uses `eps = 0.1`.

The `cost_*` outputs are sampling budgets; `eps_*` outputs are simulation errors. Check the corresponding `cvx_status_*` values.

## Related paper

Xuanqiang Zhao, Xin Wang, and Giulio Chiribella. [*Shadow Simulation of Quantum Processes*](https://doi.org/10.1103/PhysRevLett.133.120804) (2024).
