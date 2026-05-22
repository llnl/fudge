# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
# This file is required for C extensions and custom build steps.
# The main configuration is in pyproject.toml.

import glob
import os
import shutil
import subprocess
from pathlib import Path

import numpy
from setuptools import Extension, setup
from setuptools.command.build_py import build_py as _build_py

ROOT = Path(__file__).resolve().parent


class build_py(_build_py):
    """Build Python packages, and also build/stage standalone native executables."""

    def run(self):
        self._build_native_executables()
        super().run()

    def _build_native_executables(self):
        targets = [
            {
                "workdir": ROOT / "Merced",
                "source": ROOT / "Merced" / "bin" / self._exe_name("merced"),
                "dest": ROOT / "fudge" / "bin" / self._exe_name("merced"),
            },
            {
                "workdir": ROOT / "fudge" / "processing" / "deterministic" / "upscatter",
                "source": ROOT / "fudge" / "processing" / "deterministic" / "upscatter" / "bin" / self._exe_name("calcUpscatterKernel"),
                "dest": ROOT / "fudge" / "bin" / self._exe_name("calcUpscatterKernel"),
            },
        ]

        for target in targets:
            self._run_make(target["workdir"])
            self._stage_binary(target["source"], target["dest"])

    @staticmethod
    def _exe_name(name):
        return f"{name}.exe" if os.name == "nt" else name

    def _run_make(self, workdir: Path):
        subprocess.check_call(["make", "-j"], cwd=str(workdir))

    def _stage_binary(self, source: Path, dest: Path):
        if not source.exists():
            raise FileNotFoundError(f"Expected built executable not found: {source}")

        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)

        build_target = Path(self.build_lib) / dest.relative_to(ROOT)
        build_target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(dest, build_target)


setup(
    scripts=glob.glob("bin/*.py"),
    ext_modules=[
        Extension('fudge.processing.resonances._getBreitWignerSums',
            sources=['fudge/processing/resonances/getBreitWignerSums.c'], include_dirs=[numpy.get_include()], ),
        Extension('fudge.processing.resonances._getScatteringMatrices',
            sources=['fudge/processing/resonances/getScatteringMatrices.c'], include_dirs=[numpy.get_include()], ),
        Extension('fudge.processing.resonances._getCoulombWavefunctions',
            sources=['fudge/processing/resonances/getCoulombWavefunctions.c', 'fudge/processing/resonances/coulfg2.c'], include_dirs=[numpy.get_include()], ),
    ],
    install_requires=[
        'numpy',
        f'crossSectionAdjustForHeatedTarget @ {(ROOT / "crossSectionAdjustForHeatedTarget").as_uri()}',
        f'numericalFunctions @ {(ROOT / "numericalFunctions").as_uri()}',
        f'pqu @ {(ROOT / "pqu").as_uri()}',
        f'xData @ {(ROOT / "xData").as_uri()}',
        f'PoPs @ {(ROOT / "PoPs").as_uri()}',
        f'brownies @ {(ROOT / "brownies").as_uri()}',
    ],
    cmdclass={"build_py": build_py},
)
