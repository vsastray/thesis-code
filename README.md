# Numerical code for the PhD thesis: Quantum communication beyond standard quantum Shannon theory

Numerical code accompanying Xuanqiang Zhao's PhD thesis, *Quantum communication beyond standard quantum Shannon theory*.

This repository covers numerical calculations; plotting scripts are not included.

| Folder | Thesis results |
| --- | --- |
| [Chapter 3](chapter3/README.md) | Table 3.1 |
| [Chapter 4](chapter4/README.md) | Figures 4.2–4.4 |
| [Chapter 5](chapter5/README.md) | Figures 5.2–5.4 |
| [Chapter 6](chapter6/README.md) | Figure 6.2 |

## Requirements

Chapter 3 requires Python, NumPy, Paddle Quantum, and Psi4.

Chapters 4–6 require MATLAB, [CVX](https://cvxr.com/cvx/), and [QETLAB](https://qetlab.com/). Chapters 4 and 6 also require [quantinf](http://www.dr-qubit.org/download.php?file=matlab/quantinf.tar.gz) for `chanconv`. Run `cvx_setup` and add the toolboxes to MATLAB's path.

The MATLAB checks used R2025b Update 3, CVX 2.2, and QETLAB 1.0. Original dependency and solver versions are not recorded.

## Credits

Original code and numerical work are credited collectively to the authors of each chapter's associated paper, listed in the chapter READMEs. Xuanqiang Zhao assembled and adapted this collection.
