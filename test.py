from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import platform


class BuildExt(build_ext):
    def build_extensions(self):
        sys_platform = sys.platform
        machine = platform.machine()

        for ext in self.extensions:
            # --- Windows (MSVC) ---
            if self.compiler.compiler_type == "msvc":
                ext.extra_compile_args = ["/openmp"]
                ext.extra_link_args = []

            # --- Linux ---
            elif sys_platform.startswith("linux"):
                ext.extra_compile_args = ["-fopenmp"]
                ext.extra_link_args = ["-fopenmp"]

            # --- macOS (Intel or Apple Silicon) ---
            elif sys_platform == "darwin":
                # Apple Clang does not support OpenMP without libomp
                # Users must install libomp: `brew install libomp`
                ext.extra_compile_args = ["-Xpreprocessor", "-fopenmp"]
                ext.extra_link_args = ["-lomp"]

                # Special note for Apple Silicon (M1/M2/M3)
                if machine == "arm64":
                    print("Building on macOS ARM (Apple Silicon)")
                    # Optionally ensure libomp is available under correct arch
                    # You may set include_dirs and library_dirs manually if needed
                    # ext.include_dirs = ['/opt/homebrew/include']
                    # ext.library_dirs = ['/opt/homebrew/lib']

            else:
                print(f"Warning: Unsupported platform: {sys_platform}")

        super().build_extensions()


if __name__ == "__main__":
    BuildExt()
