# OPenRCWA

[![CI](https://github.com/sirokujira/openrcwa/actions/workflows/ci.yml/badge.svg)](https://github.com/sirokujira/openrcwa/actions/workflows/ci.yml)

A free and open source RCWA solver, built alongside an OpenFDTD-derived FDTD solver.

## Continuous Integration

`.github/workflows/ci.yml` builds the project (CPU-only; CUDA/MKL/MPI disabled)
on Ubuntu and runs a lightweight smoke test on every push and pull request:

- Configures and builds the FDTD solver (`orcwa`, `orcwa_post`), the RCWA
  solver test executables, and the C data-creation library + sample.
- Generates a sample input file with both the C (`sample_rcwa_grating`) and
  Python (`python/datalib/sample_rcwa_grating.py`) data-creation libraries and
  checks the Python output against the committed file.
- Runs the solver on a reduced smoke-test input (`ci/smoke_test.orcwa`) and
  verifies it reaches a normal end.

To reproduce the CI build locally:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF
cmake --build build -j"$(nproc)"
```
