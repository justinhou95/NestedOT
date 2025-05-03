from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
import os
import platform
import subprocess

__version__ = "0.0.1"

compile_args = ["-std=c++14"]
link_args = []

# Detect platform specifics
system = platform.system()

# Enable OpenMP if available
if system == "Darwin":
    compile_args += ["-Xpreprocessor", "-fopenmp"]
    try:
        # Try to get the libomp prefix via brew
        prefix = (
            subprocess.check_output(["brew", "--prefix", "libomp"]).decode().strip()
        )
        libomp_lib_dir = os.path.join(prefix, "lib")
    except Exception:
        # Fallback if brew command fails
        libomp_lib_dir = "/usr/local/opt/libomp/lib"
    link_args += [f"-L{libomp_lib_dir}", "-lomp"]

    openmp_include_dir = "/usr/local/opt/libomp/include"  # Adjust if on Apple Silicon, e.g., "/opt/homebrew/opt/libomp/include"
    extra_include_dirs = [openmp_include_dir]
elif system == "Linux":
    compile_args += ["-fopenmp"]
    link_args += ["-fopenmp"]
    extra_include_dirs = []
elif system == "Windows":
    compile_args += ["/openmp"]
    extra_include_dirs = []
else:
    extra_include_dirs = []

# Set up include directories
include_dirs = ["src/include"]
include_dirs += extra_include_dirs

# Use environment variable to get Eigen include directory if provided.
eigen_dir = os.environ.get("EIGEN3_INCLUDE_DIR")
if eigen_dir:
    include_dirs.append(eigen_dir)
else:
    # Fallback: check common eigen installation directories.
    common_eigen_dirs = ["/usr/include/eigen3", "/usr/local/include/eigen3"]
    for d in common_eigen_dirs:
        if os.path.exists(d):
            include_dirs.append(d)
            break
    else:
        print(
            "Warning: Eigen include directory not found. Please set EIGEN3_INCLUDE_DIR."
        )

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
        extra_compile_args=compile_args,
        extra_link_args=link_args,
        cxx_std=14,
    ),
]

setup(
    name="pnot",
    version=__version__,
    description="Adapted Wasserstein distance solver in C++/pybind11",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
)
