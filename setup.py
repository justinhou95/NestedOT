# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages
import os
import platform

__version__ = "0.0.1"

if platform.system() == "Darwin":
    os.environ["CC"] = "/opt/homebrew/opt/llvm/bin/clang"
    os.environ["CXX"] = "/opt/homebrew/opt/llvm/bin/clang++"

ext_modules = [
    Pybind11Extension(
        "_wrapper",
        [
            "src/wrapper.cpp",
            "src/solver.cpp",
            "src/emd_wrap.cpp",
            "src/utils.cpp",
            "src/printer.cpp",
        ],
        # Example: passing in the version to the compiled code
        define_macros=[("VERSION_INFO", __version__)],
        include_dirs=["extern/eigen", "./src/include"],
        extra_compile_args=["-fopenmp"],
        extra_link_args=["-fopenmp"],
    ),
]

setup(
    name="pnot",
    version=__version__,
    description="Nested Optimal Transport",
    long_description="",
    ext_modules=ext_modules,
    packages=find_packages(),
    extras_require={"test": "pytest"},
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
)
