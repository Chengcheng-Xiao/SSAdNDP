---
anchor: visualization
title: Visualization
order: 6
---

## Visualization

If `T` is entered on line 3 of the input file, the program creates <code>vis_gs_<em>n</em>c-bonds.out</code> files (one per bond type, general search) or `vis_ud_bonds.out` (user-directed). For spin-polarized calculations, separate files are created for each spin channel.

Use `visual.exe` to convert these to Gaussian cube format:

```sh
./path-to-Visualization/visual.exe vis_gs_3c-bonds.out
```

Three optional arguments specify grid resolution (number of integration points per cell direction; default: 50 in each):

```sh
./path-to-Visualization/visual.exe vis_gs_3c-bonds.out 50 50 50
```

By default, bonds are represented in a supercell containing the central and all neighboring cells (since bonds may span unit cells). To restrict to the central cell only (useful for a finite molecule in a large cell), append `F`:

```sh
./path-to-Visualization/visual.exe vis_gs_3c-bonds.out 50 50 50 F
```

### VESTA Workflow

<div class="vesta-steps" markdown="1">
1. Open `lattice_vec.cube` in [VESTA](https://jp-minerals.org/vesta/en/). Uncheck *Show sections*, check *Show isosurfaces*.
2. Import bonds: **Edit → Edit Data → Phase → Import**, select a `*_bond_*.cube` file.
3. Select the first entry and set the three numbers under *This layer* to `0.5 0.5 0.5` to align the supercell origin correctly.
</div>
