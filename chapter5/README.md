# Chapter 5 — Figures 5.2–5.4

Requires MATLAB, CVX, and QETLAB; see [setup](../README.md#requirements).

## Run

Set this folder as MATLAB's current folder and run the relevant script.

| Figure | Script | Calculation |
| --- | --- | --- |
| 5.2 | [info_recovery.m](info_recovery.m) | Recovery costs for three channels at noise levels `0:0.1:0.5` |
| 5.3 | [data_generating.m](data_generating.m) | Repeated sampling of observable estimates |
| 5.4 | [entry_extraction.m](entry_extraction.m) | TC and QPD costs; select `case_id` from 1 to 7 |

Channel labels are `DE` (depolarizing), `AD` (amplitude damping), and `DP` (dephasing).

Figure 5.3 writes `sample_*.csv` in the current folder, overwriting existing files. It preserves the original observable convention `Y = [0 1i; -1i 0]`.

Figure 5.4 uses a reconstructed driver: the original fixed-example driver was unavailable. Its seven patterns match the corrected thesis figure, and its costs agree numerically with the saved plotted values.

## Related paper

Xuanqiang Zhao, Lei Zhang, Benchi Zhao, and Xin Wang. [*Power of quantum measurement in simulating unphysical operations*](https://doi.org/10.1103/PhysRevResearch.7.013334) (2025).
