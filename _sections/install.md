---
anchor: install
title: Installation
order: 2
---

## Installation

The package has three components: a DFT interface + projection executable, the SSAdNDP calculation executable, and a visualization executable. The code is written in Fortran 90 and requires `BLAS95` and `LAPACK95` as implemented in Intel MKL, plus OpenMP.

### DFT Interface — VASP

Compatible with VASP v5.4.4.

1. Copy `interfaces/VASP_interface/projection_output.F` into the VASP source directory (`./scr`).
2. Apply the patch from the VASP root directory:

   ```sh
   patch -p0 < proj_VASP_5.4.patch
   ```

3. Compile VASP as usual.

### DFT Interface — CRYSTAL

No modification needed.

### SSAdNDP Executable

```sh
# 1. Edit SSAdNDP/Makefile.template — set MKLROOT
./mkmf -t Makefile.template -p ssadndp.exe
make
```

### Projection Executable

```sh
# Edit projection/Makefile — set MKLROOT
make
```

### Visualization Executable

```sh
ifort visual.f90 -o visual.exe
```
