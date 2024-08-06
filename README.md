# climlab-rrtmg

[![Build and test](https://github.com/climlab/climlab-rrtmg/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/climlab/climlab-rrtmg/actions/workflows/build-and-test.yml)

Brian Rose, University at Albany

## About

This is a stand-alone Python wrapper for the [RRTMG](http://rtweb.aer.com/rrtm_frame.html) radiation modules.

The primary use-case is to drive the RRTMG radiation processes in [climlab](https://climlab.readthedocs.io/),
but it can be used as a stand-alone radiation model if you are familiar with the
RRTMG Fortran interface. This is a lightweight wrapper that emulates the Fortran
interface as closely as possible.

Currently we are wrapping RRTMG_LW v4.85 and RRTMG_SW v4.0. The Fortran source code
is included in this repository. The latest versions 5.0 of the RRTMG source code
are available on GitHub [here](https://github.com/AER-RC)

This wrapper includes a modification to RRTMG_LW to allow reporting of OLR
components in spectral bands, as illustrated using climlab
[here](https://climlab.readthedocs.io/en/latest/courseware/Spectral_OLR_with_RRTMG.html).
This modification is strictly diagnostic does not change any other behavior of RRTMG_LW.

## Installation

Pre-built binaries for many platforms are available from [conda-forge](https://conda-forge.org).

To install in the current environment:
```
conda install climlab-rrtmg --channel conda-forge
```
or create a self-contained environment:
```
conda create --name my_env python=3.10 climlab-rrtmg --channel conda-forge
conda activate my_env
```

See below for instructions on how to build from source.

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

### Build environment

Here are instructions to create a build environment (including Fortran compiler) with conda/mamba

Starting from the root of the `climlab-rrtmg` repo *(example for Apple M1 machine, see `./ci/` for other environment files)*:
```
mamba env create --file ./ci/requirements-macos-arm64.yml
conda activate rrtmg_build_env
```

Or, to specify the Python version, you can do
```
mamba create --name rrtmg_build_env python=3.11 --channel conda-forge
mamba env update --file ./ci/requirements-macos-arm64.yml
conda activate rrtmg_build_env
```

### Building and installing into the Python environment

From the root of the repository, do this:
```
python -m pip install . --no-deps -vv
```

### Running tests

To run tests, do this from any directory other than the climlab-emanuel-convection repo:
```
pytest -v --pyargs climlab_sbm_convection
```

## Version history

Version 0.3 updates the build infrastructure.

Version 0.2 is the first public release (April 2022).
The Python wrapper code has been extracted from
[climlab v0.7.13](https://github.com/brian-rose/climlab/releases/tag/v0.7.13).
