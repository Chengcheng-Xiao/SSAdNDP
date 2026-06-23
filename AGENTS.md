# SSAdNDP — Agent Guide

## Branches

- **`master`** — Fortran 90 source code (SSAdNDP/, projection/, visualization/, interfaces/)
- **`site`** (default) — GitHub Pages site (Jekyll project: `_config.yml`, `Gemfile`, `_layouts/`, `assets/`, `index.html`)

Always work on `master` for code changes.

## Build (master branch)

All executables require Intel Fortran (`ifort`) + Intel MKL. Set `MKLROOT` in the relevant makefile first.

```sh
# SSAdNDP executable
cd SSAdNDP && ./mkmf -t Makefile.template -p ssadndp.exe && make

# Projection executable
cd projection && make

# Visualization executable
cd visualization && ifort visual.f90 -o visual.exe
```

No tests, no linter, no formatter, no typechecker exist. No CI.

## Workflow

1. Run VASP (v5.4.4) with `LNBO=.T.`, `ISYM=0`, default `NPAR` → produces `wavefunction.dat`
2. Project: `projection.exe basis.inp wavefunction.dat NBO.out`
3. SSAdNDP general search: `ssadndp.exe SSAdNDP.inp`
4. SSAdNDP user-directed: `ssadndp.exe SSAdNDP-UD.inp`
5. Visualize: `visual.exe vis_ud_bonds.out`

Run with `ulimit -s unlimited` for large systems.

## Key source files

| File | Purpose |
|------|---------|
| `SSAdNDP/nbo_main.f90` | Entry point; reads NBO.out, calls SSAdNDP module |
| `SSAdNDP/ssadndp_module.f90` | General + user-directed search logic |
| `projection/projection_main.f90` | Wavefunction → AO projection |
| `visualization/visual.f90` | Bond output → Gaussian cube format |

## Input files

- `SSAdNDP.inp` — general search: bond types (n-centers), ON threshold, distance threshold
- `SSAdNDP-UD.inp` — user-directed: explicit atom+cell lists per fragment
- `CONTCAR` — read for lattice vectors and atomic coordinates (must use Direct coordinates)
- `NBO.out` / `NBO_mat.out` — projection output read by SSAdNDP

## Checkpoint / resume

```sh
# Write checkpoint
ssadndp.exe SSAdNDP.inp chkfilename
# Read checkpoint without overwriting
ssadndp.exe SSAdNDP.inp chkfilename F
```

Residual density from previous run should match initial density of resumed run.

## Spin-polarized

Input files for ISPIN=2 have separate sections for spin 1 and spin 2. Visual output gets `_spin-1_` / `_spin-2_` suffixes.

## Visualization

`visual.exe` converts `vis_*_bonds.out` to `.cube` files. Default grid: 50×50×50 points per cell direction (overridable via CLI args). Append `F` to restrict to central cell only.
