# climlab-rrtmg

![GitHub Workflow Status (branch)](https://img.shields.io/github/workflow/status/brian-rose/climlab-rrtmg/build-and-test/main?logo=github&style=for-the-badge)

Brian Rose, University at Albany

## About

This is a stand-alone Python wrapper for the RRTMG radiation modules.

The primary use-case is to drive the RRTMG radiation processes in [climlab](https://climlab.readthedocs.io/),
but it can be used as a stand-alone radiation model if you are familiar with the
RRTMG Fortran interface. This is a lightweight wrapper that emulates the Fortran
interface as closely as possible.

## Installation

Pre-built binaries for many platforms will be available to install from conda-forge. Stay tuned.

For now, please build from source (instructions below).

## Example usage

You can import the modules into a Python session with
```
from climlab_rrtmg import rrtmg_lw, rrtmg_sw
```

The main RRTMG drivers are exposed through
```
rrtmg_lw.climlab_rrtmg_lw()
```
and
```
rrtmg_sw.climlab_rrtmg_sw()
```

Please see the directory `climlab_rrtmg/tests/` directory in this repository
for working examples that set up all the necessary input arrays and call the drivers.

## Building from source

Here are instructions to create a build environment (including Fortran compiler)
with conda/mamba and build using f2py.
It should be possible to build using other Fortran compilers, but I haven't tested this.

To build *(example for Apple M1 machine, see `./ci/` for other environment files)*:
```
mamba create --name rrtmg_build_env python=3.10 --channel conda-forge
mamba env update --file ./ci/requirements-macos-arm64.yml
conda activate rrtmg_build_env
python -m pip install . --no-deps -vv
```

To run tests, do this from any directory other than the climlab_rrtmg repo:
```
pytest -v --pyargs climlab_rrtmg
```
