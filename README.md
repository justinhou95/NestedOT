# PNOT: Python Nested Optimal Transport 🪆

This library implements very fast C++ and Python solver of nested (adapted) optimal transport problem. In particular, it calculates the adapted Wasserstein distance fast and accurately. We provide both C++ and Python implementation, and a wrapper to use the fast C++ solver with Python. This solver is very easy to use and all you need to do is feeding two paths distribution into the solver. The solver will do all the adapted empirical measures, quantization, and nested computation for you automatically and swiftly.

## Installation 📦

### Preparation for Mac User (No need for Win/Linux)

For Mac user, make sure you install the Apple’s default clang compiler (Xcode):
```bash
$ xcode-select --install
```
Then we need to install the official LLVM clang/clang++ compilers (because Apple builds its version of clang without OpenMP support) with Homebrew:

```bash
$ brew install llvm
$ brew install libomp
```

### Installation
There are few options for installation:

Install the stable version via pip:
```bash
$ pip install pnot
``` 
Install the latest github version:
```bash
$ pip install git+https://github.com/justinhou95/NestedOT.git
``` 
Clone the github repo for development to access to tests, tutorials and scripts and install in [develop mode](https://setuptools.pypa.io/en/latest/userguide/development_mode.html):
```bash
$ git clone https://github.com/justinhou95/NestedOT.git
```
```bash
$ cd NestedOT
$ pip install -e .
``` 

## Notebooks
- [demo.ipynb](https://github.com/justinhou95/NestedOT/blob/main/notebooks/demo.ipynb) shows quickstart using PNOT to compute nested distance
- [solver_explain.ipynb](https://github.com/justinhou95/NestedOT/blob/main/notebooks/solver_explain.ipynb) explains how the conditional distributions are estimated and how nested computing works.
- [exemple_of_use.ipynb](https://github.com/justinhou95/NestedOT/blob/main/notebooks/exemple_of_use.ipynb) shows how to use it in practice.


## Reference
* [Fast Transport](https://github.com/nbonneel/network_simplex/tree/master)
* [Python Optimal Transport](https://github.com/PythonOT/POT)
* [Entropic adapted Wasserstein distance on Gaussians](https://arxiv.org/abs/2412.18794)
