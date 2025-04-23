# NestedOT  
**Nested Optimal Transport** (a.k.a. Adapted Wasserstein)  
Fast, exact C++ solver for discrete nested OT via backward induction—unmatched performance.  

---

## Why NestedOT?  
The only other publicly available implementations of exact discrete adapted OT (e.g. [stephaneckstein/aotnumerics](https://github.com/stephaneckstein/aotnumerics)) become utterly infeasible beyond very small sample sizes. **NestedOT** closes that gap, delivering:

- **→ Over 1000× speed-up** on non-Markovian problems  **TO BE COMPLETED: Markovian benchmarks**  
- **→ Scales to thousands of samples** in seconds  
- **→ C++ core with Python bindings** for seamless integration  

---

## Performance Comparison: Full-History OT  

We compare **NestedOT’s** C++ `nested_ot` solver against the pure-Python dynamic-programming reference (`solve_dynamic`) for the **non-Markovian** case.

- **Problem setup**  
  - Sample sizes: **100, 200, 300, 500, 1000**  
  - **5 runs** each  
  - Matrices  
    \[
      L = \begin{pmatrix}
        1 & 0 & 0 & 0\\
        2 & 2 & 0 & 0\\
        1 & 1 & 3 & 0\\
        2 & 2 & 1 & 2
      \end{pmatrix},\quad
      M = \begin{pmatrix}
        1 & 0 & 0 & 0\\
        2 & 1 & 0 & 0\\
        3 & 2 & 1 & 0\\
        4 & 3 & 2 & 1
      \end{pmatrix}
    \]
  - Grid size \(\delta = 0.2\), cost exponent \(p=2\)  
  - **Backward induction** on the full-history DAG  

| Sample \(n\) | C++ time (mean ± std) [s] | Py DP time (mean ± std) [s] | Speed-up (×)    |
|-------------:|:-------------------------:|:--------------------------:|:---------------:|
|       **100**  |     0.004 ± 0.000         |     3.753 ± 0.124          | **1067.7×**     |
|       **200**  |     0.011 ± 0.000         |    14.746 ± 0.397          | **1394.4×**     |
|       **300**  |     0.019 ± 0.001         |    33.338 ± 0.782          | **1737.3×**     |
|       **500**  |     0.045 ± 0.004         |    92.299 ± 1.687          | **2070.6×**     |
|      **1000**  |     0.165 ± 0.002         |   405.226 ± 6.960          | **2460.7×**     |

![Timing vs. Sample Size for Full-History OT](./full_history_timing.png)

> **Over three orders of magnitude** faster than any existing public code—and the gap **widens** with larger samples.

---

## Getting Started  
**TO BE COMPLETED:** installation instructions, example usage, etc.
