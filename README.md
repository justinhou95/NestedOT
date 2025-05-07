# PNOT: Python Nested Optimal Transport 🪆

This library implements a very fast C++ and Python solver for the nested (adapted) optimal transport problem. In particular, it calculates the adapted Wasserstein distance quickly and accurately. We provide both C++ and Python implementations, and a wrapper to use the fast C++ solver from Python. Just feed two path distributions into the solver—the rest (empirical measures, quantization, nested computation) happens automatically and swiftly.

## Installation 📦

### Preparation for macOS Users
Make sure you have Apple’s Xcode command-line tools installed:
```bash
$ xcode-select --install
```
Install LLVM and OpenMP support via Homebrew:
```bash
$ brew install llvm libomp
```

### Installation

- **Stable release** via PyPI:
  ```bash
  $ pip install pnot
  ```
- **Latest GitHub version**:
  ```bash
  $ pip install git+https://github.com/justinhou95/NestedOT.git
  ```
- **Developer mode** (clone and install editable):
  ```bash
  $ git clone https://github.com/justinhou95/NestedOT.git
  $ cd NestedOT
  $ pip install -e .
  ```

## Notebooks

- [demo.ipynb](https://github.com/justinhou95/NestedOT/blob/main/notebooks/demo.ipynb) — Quickstart and basic usage
- [solver_explain.ipynb](https://github.com/justinhou95/NestedOT/blob/main/notebooks/solver_explain.ipynb) — How conditional distributions are estimated and nested computations performed
- [example_of_use.ipynb](https://github.com/justinhou95/NestedOT/blob/main/notebooks/exemple_of_use.ipynb) — Approach similar to that described in Backhoff et al. 2021 for estimating adapted Wasserstein distance with continuous measures

## Performance Comparison: Non-Markovian

We compare **PNOT’s** C++ `nested_ot` solver against the only publicly available alternative—`solve_dynamic` from AOTNumerics (Eckstein & Pammer 2023)—for the **non-Markovian** case.

- **Sample sizes**: 100, 200, 300, 500, 1000 (5 runs each)  
- **Matrices**: L and M (shown below) are lower-triangular Cholesky factors used to generate zero-mean Gaussian path distributions:

$$
L = \begin{pmatrix} 1 & 0 & 0 & 0\\ 2 & 2 & 0 & 0\\ 1 & 1 & 3 & 0\\ 2 & 2 & 1 & 2 \end{pmatrix}, \quad M = \begin{pmatrix} 1 & 0 & 0 & 0\\ 2 & 1 & 0 & 0\\ 3 & 2 & 1 & 0\\ 4 & 3 & 2 & 1 \end{pmatrix} 
$$

- **Grid size**: \(\delta=0.2\), **cost exponent**: \(p=2\)

| Sample \(n\) | C++ time (mean ± std) [s] | AOTNumerics `solve_dynamic` [s] | Speed‑up (×)  |
|-------------:|:-------------------------:|:-------------------------------:|:-------------:|
| 100          | 0.040 ± 0.011             | 1.866 ± 0.061                   | 46.5×         |
| 200          | 0.024 ± 0.010             | 7.490 ± 0.194                   | 317.5×        |
| 300          | 0.015 ± 0.004             | 17.061 ± 0.457                  | 1124.9×       |
| 500          | 0.019 ± 0.000             | 48.557 ± 0.897                  | 2543.1×       |
| 1000         | 0.067 ± 0.014             | 221.924 ± 4.303                 | 3292.6×       |

![Timing vs. Sample Size for Full-History OT](./full_history_timing.png)

> **Over three orders of magnitude** faster than the only other public implementation—and the gap widens with larger samples.

## Performance Comparison: Markovian

We compare **PNOT’s** C++ `nested_ot` solver against AOTNumerics’ `solve_dynamic` for the **Markovian** case.

- **Time steps**: \(T=10\)  
- **Sample sizes**: 100, 200, 300, 500, 1000 (5 runs each)
- **Grid size**: \(\delta=0.1\), **cost exponent**: \(p=2\)  
- **Path generation**: integrated-process (variance 0.25) & AR(1) (\(\phi=0.8,\sigma=1.0\))

| Sample \(n\) | C++ time (mean ± std) [s] | AOTNumerics `solve_dynamic` [s] | Speed‑up (×)  |
|-------------:|:-------------------------:|:-------------------------------:|:-------------:|
| 200          | 0.007 ± 0.000             | 8.531 ± 0.187                   | 1224.3×       |
| 600          | 0.022 ± 0.001             | 46.357 ± 0.982                  | 2120.2×       |
| 1000         | 0.036 ± 0.000             | 95.701 ± 0.924                  | 2637.6×       |
| 1500         | 0.052 ± 0.001             | 157.432 ± 2.981                 | 3002.5×       |
| 2000         | 0.069 ± 0.001             | 221.344 ± 4.013                 | 3215.8×       |

![Timing vs. Sample Size for Markovian OT](./Markov_timing.png)

## Reference

- Eckstein & Pammer, *Computational methods for adapted optimal transport*, 2023, arXiv:2203.05005 ([PDF](https://arxiv.org/abs/2203.05005)) — only other public solver (`solve_dynamic`)  
- Backhoff et al., *Estimating processes in adapted Wasserstein distance*, 2021, arXiv:2002.07261 ([PDF](https://arxiv.org/abs/2002.07261)) — continuous-measure discretization strategy  
- [Fast Transport (Network Simplex)](https://github.com/nbonneel/network_simplex/tree/master)  
- [Python Optimal Transport (POT)](https://github.com/PythonOT/POT)  
- [Entropic Adapted Wasserstein on Gaussians](https://arxiv.org/abs/2412.18794)  

