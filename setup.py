import os
import re
import sys
import subprocess
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        cmake_executable = os.path.join(sys.prefix, "bin", "cmake")
        if not os.path.exists(cmake_executable):
            cmake_executable = "cmake"
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name))
        )

        cfg = "Release"

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPython3_EXECUTABLE={sys.executable}",
            f"-DCMAKE_PREFIX_PATH={sys.prefix}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
        ]

        build_args = ["--config", cfg]

        build_temp = self.build_temp
        os.makedirs(build_temp, exist_ok=True)

        subprocess.check_call(
            [cmake_executable, ext.sourcedir] + cmake_args,
            cwd=build_temp,
        )

        subprocess.check_call(
            [cmake_executable, "--build", "."] + build_args,
            cwd=build_temp,
        )


setup(
    name="qdex",
    version="1.0.0",
    description="QDEX: Quantum Dot Excitations & Dynamics",
    author="Ivan Infante et al.",
    packages=find_packages(),
    ext_modules=[CMakeExtension("libint_cpp", sourcedir="libint")],
    cmdclass={"build_ext": CMakeBuild},
    entry_points={
        "console_scripts": [
            "qdex=qdex.cli:main",
            "minibse=qdex.cli:main",
        ],
    },
    zip_safe=False,
)
