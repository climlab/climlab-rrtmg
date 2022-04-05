# climlab_rrtmg

Just trying out a more federated approach to climlab submodules

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
