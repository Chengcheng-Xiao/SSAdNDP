---
anchor: userdirected
title: User-Directed Search
order: 5
---

## User-Directed Search

For higher *n*, the combinatorial space becomes intractable. The user-directed mode lets you probe specific sets of atomic centers for bonding. This may introduce bias relative to general search, but enables chemically meaningful results on fragments of interest. It is usually run after general search, reading the residual density matrix from a checkpoint file.

<div class="note" markdown="1">
**Note:** Unlike general search (which depletes density after all bonds of a given *n* are found), user-directed search updates the density matrix after every single bond. For adjacent symmetrically equivalent bonds, resulting ONs may differ slightly; the averaged value is typically accepted.
</div>

### Running

```sh
./path-to-SSAdNDP/ssadndp.exe SSAdNDP-UD.inp chkpoint_from_gen F
```

### Input File — Non-Polarized

<div class="input-block" markdown="1">
<div class="input-block-label">SSAdNDP-UD.inp</div>

```
USER-DIRECTED SEARCH INPUT FILE
Prepare files for visualization: T/F?
F
Number of fragments to search on:
2
Number of centers on fragment 1:
6
Atom numbers, cell a,b,c:
1, 0 0 0
2, 0 0 0
7, 0 -1 0
5, 0 -1 0
6, -1 0 0
4, -1 0 0
Number of bonds on fragment 1:
1
Number of centers on fragment 2:
7
Atom numbers, cell a,b,c:
1, 0 0 0
2, 0 0 0
3, 0 0 0
4, 0 0 0
5, 0 0 0
6, 0 0 0
7, 0 0 0
Number of bonds on fragment 2:
1
```

</div>

Cell indices are relative to the central cell (0,0,0). An atom listed as `2, 1 0 -1` refers to the image of atom 2 obtained by translation of (+**a** − **c**).

### Input File — Spin Polarized

<div class="input-block" markdown="1">
<div class="input-block-label">SSAdNDP-UD.inp (spin-polarized)</div>

```
USER-DIRECTED SEARCH INPUT FILE
Prepare files for visualization: T/F?
T
Number of fragments to search on (spin 1):
1
Number of centers on fragment 1:
6
Atom numbers, cell a,b,c:
1, 0 0 0
2, 0 0 0
3, 0 0 0
4, 0 0 0
5, 0 0 0
6, 0 0 0
Number of bonds on fragment 1:
3
Number of fragments to search on (spin 2):
1
Number of centers on fragment 1:
6
Atom numbers, cell a,b,c:
1, 0 0 0
2, 0 0 0
3, 0 0 0
4, 0 0 0
5, 0 0 0
6, 0 0 0
Number of bonds on fragment 1:
2
```

</div>
