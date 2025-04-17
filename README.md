# PNOT Python Nested Optimal Transport 🚀

This library implements very fast C++ and Python solver of nested (adapted) optimal transport problem. In particular, it calculate the adapted Wasserstein distance fast and accurately. We provide both C++ and Python implementation, and a wrapper to use the fast C++ solver with Python. This solver is very easy to use and all you need to do if feeding two paths distribution into the solver function and the solver will do all the adapted empirical measures, quantization, and nested computation for you.

[documentation](https://justinhou95.github.io/NeuralHedge/)

## Installation 📦

Install the latest stable release:

```bash
$ pip install pnot==1.0.0
``` 

Install the latest github version:

```bash
$ pip install git+https://github.com/justinhou95/NestedOT.git
``` 

Clone the github repo for development to access to tests, tutorials and scripts.
```bash
$ git clone https://github.com/justinhou95/NestedOT.git
```
and install in [develop mode](https://setuptools.pypa.io/en/latest/userguide/development_mode.html)
```bash
$ cd NestedOT
$ pip install -e .
``` 

## Tutorials 📚
To help you to understand how easy to train your models with PNOT, we also provide tutorials in the notebooks folder.

## Dealing with issues 🛠️

If you are experiencing any issues while running the code or request new features/models to be implemented please [open an issue on github](https://github.com/justinhou95/NestedOT/issues).


## Contributing 🚀
You want to contribute to this library, that will be cool!
