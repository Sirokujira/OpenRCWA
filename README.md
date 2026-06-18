# OpenRCWA

[![CI](https://github.com/sirokujira/openrcwa/actions/workflows/ci.yml/badge.svg)](https://github.com/sirokujira/openrcwa/actions/workflows/ci.yml)

A free and open source RCWA solver, built alongside an OpenFDTD-derived FDTD solver.

## Solvers

The project shares a single text input format (`.orcwa`) across two solvers:

- **`orcwa`** — the OpenFDTD-derived FDTD time-domain solver.
- **`orcwa_rcwa`** — a frequency-domain RCWA solver for layered periodic
  structures (diffraction gratings, thin-film stacks, metasurfaces).

### RCWA solver (`orcwa_rcwa`)

`rcwa/RCWAInput.{h,cpp}` parses the same `.orcwa` file used by the FDTD
solver and interprets it as a layered periodic structure:

- `xmesh`/`ymesh` extents define the in-plane periods.
- `zmesh` and `geometry` boxes are sliced into homogeneous/patterned layers
  along z. The top (largest z) and bottom (smallest z) slabs are treated as
  the semi-infinite incidence and transmission media.
- `material` entries set each region's permittivity; `planewave` sets the
  incidence; `frequency1` is swept and converted to wavelengths.
- `rcwaorder = nHx [nHy]` sets the number of Fourier harmonics
  (total `2*nH+1`); the default is `nHx = 10`, `nHy = 0` (1D gratings).

```sh
orcwa_rcwa [-o <out.csv>] <datafile.orcwa>
```

It writes a `lambda, R, T, R+T` spectrum (reflectance / transmittance summed
over all diffraction orders). The implementation is validated against analytic
references: an air→glass interface (Fresnel `R = 0.04`), a Fabry–Pérot slab
(energy-conserving interference extrema), and a dielectric grating (energy
conservation with the expected Rayleigh anomaly).

> Note: the current RCWA driver targets normal-incidence 1D gratings
> (`nHy = 0`); oblique incidence and full 2D periodicity are future work.

## Continuous Integration

`.github/workflows/ci.yml` builds the project (CPU-only; CUDA/MKL/MPI disabled)
on Ubuntu and runs a lightweight smoke test on every push and pull request:

- Configures and builds the FDTD solver (`orcwa`, `orcwa_post`), the RCWA
  solver test executables, and the C data-creation library + sample.
- Generates a sample input file with both the C (`sample_rcwa_grating`) and
  Python (`python/datalib/sample_rcwa_grating.py`) data-creation libraries and
  checks the Python output against the committed file.
- Runs the FDTD solver on a reduced smoke-test input (`ci/smoke_test.orcwa`)
  and verifies it reaches a normal end.
- Runs the RCWA solver on an air→glass Fresnel case (`ci/smoke_rcwa.orcwa`)
  and checks `R ≈ 0.04` and `R + T ≈ 1` (energy conservation).

To reproduce the CI build locally:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DWITH_CUDA=OFF
cmake --build build -j"$(nproc)"
```
