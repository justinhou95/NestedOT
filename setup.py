from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages
import os
import platform
import sys

__version__ = "0.0.1"

# Detect platform and arch
system = platform.system()
machine = platform.machine()

# Start with defaults
extra_compile_args = []
extra_link_args = []

if system == "Windows":
    # MSVC compiler (e.g., Visual Studio)
    extra_compile_args = ["/openmp"]
    extra_link_args = []

elif system == "Linux":
    extra_compile_args = ["-fopenmp"]
    extra_link_args = ["-fopenmp"]

elif system == "Darwin":
    extra_compile_args = ["-Xpreprocessor", "-fopenmp"]
    extra_link_args = ["-fopenmp", "-lomp"]

    # macOS needs extra handling for Apple Clang
    # Requires: brew install llvm libomp
    # Use Homebrew's LLVM for proper OpenMP support

    if machine == "arm64":
        print("Detected macOS on Apple Silicon (M1/M2/M3)")
        llvm_root = "/opt/homebrew/opt/llvm"
    elif machine == "x86_64":
        print("Detected macOS on Intel")
        llvm_root = "/usr/local/opt/llvm"
    else:
        raise RuntimeError(f"Unsupported machine: {machine}")

    os.environ["CC"] = f"{llvm_root}/bin/clang"
    os.environ["CXX"] = f"{llvm_root}/bin/clang++"

else:
    raise RuntimeError(f"Unsupported platform: {system}")

include_dirs = ["extern/eigen", "./src/include"]

ext_modules = [
    Pybind11Extension(
        "_wrapper",
        sources=[
            "src/wrapper.cpp",
            "src/solver.cpp",
            "src/emd_wrap.cpp",
            "src/utils.cpp",
            "src/printer.cpp",
        ],
        define_macros=[("VERSION_INFO", __version__)],
        include_dirs=include_dirs,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    )
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
