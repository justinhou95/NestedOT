# PNOT: Python Nested Optimal Transport 🪆

This library implements very fast C++ and Python solver for the nested (adapted) optimal transport problem. In particular, it calculates the adapted Wasserstein distance quickly and accurately. We provide both C++ and Python implementations, and a wrapper to use the fast C++ solver from Python. It is very easy to use. Just feed two path samples into the solver—the rest (empirical measures, quantization, nested computation) happens automatically and swiftly.

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
- [example_of_use.ipynb](https://github.com/justinhou95/NestedOT/blob/main/notebooks/example_of_use.ipynb) — Approach similar to that described in Backhoff et al. 2021 for estimating adapted Wasserstein distance with continuous measures

## Performance Comparison:
We compare **PNOT’s** C++ `nested_ot` solver against the only publicly available alternative—`solve_dynamic` from AOTNumerics (Eckstein & Pammer 2023). For both markovian and non-markovian solver, we obtain more than 3000× speed improvement and the gap widens with larger samples. 

![Timing vs. Sample Size for Full-History OT](./assets/Markovian.png)

## Convergence Analysis:
We test the consistency of our solver as number of samples increases. By choosing the grid size $\Delta_N = N^{-\frac{1}{dT}}$, the adapted 2-Wasserstein distance between adapted empirical measures converge to the theoretical adapted 2-Wasserstein distance between underlying measures, namely $\mathcal{A}\mathcal{W}_2(\hat{\mu}^N, \hat{\nu}^N) \to \mathcal{A}\mathcal{W}_2(\mu, \nu)$. This numerically confirm Theorem 2.7 in [1].

![Convergence](./assets/ConvergenceMarkovian.jpg)



## Reference

- [1] Acciaio, B., & Hou, S. (2024). Convergence of adapted empirical measures on R d. The Annals of Applied Probability, 34(5), 4799-4835. ([PDF](https://arxiv.org/pdf/2211.10162))
- [2] Eckstein, S., & Pammer, G. (2024). Computational methods for adapted optimal transport. The Annals of Applied Probability, 34(1A), 675-713. ([PDF](https://arxiv.org/abs/2203.05005)) — only other public solver (`solve_dynamic`) 
- [3] Backhoff, J., Bartl, D., Beiglböck, M., & Wiesel, J. (2022). Estimating processes in adapted Wasserstein distance. The Annals of Applied Probability, 32(1), 529-550. ([PDF](https://arxiv.org/abs/2002.07261)) — continuous-measure discretization strategy  
- [Fast Transport (Network Simplex)](https://github.com/nbonneel/network_simplex/tree/master)  
- [Python Optimal Transport (POT)](https://github.com/PythonOT/POT)  
- [Entropic Adapted Wasserstein on Gaussians](https://arxiv.org/abs/2412.18794)  

