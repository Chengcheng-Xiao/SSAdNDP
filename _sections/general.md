---
anchor: general
title: General Search
order: 4
---

## General Search

The general search algorithm tests all possible combinations of *n* atoms for *n*c-2e bonds. The user specifies bond types (by number of centers) and occupation number (ON) thresholds. A bond is accepted if its ON exceeds the threshold (typically close to 2.00 \|e\|). After all combinations of *n* atoms are tested, accepted bonds are depleted from the density matrix and the search continues for (*n*+1)c-2e bonds.

### Running

```sh
./path-to-SSAdNDP/ssadndp.exe SSAdNDP.inp
```

Two optional command-line arguments: a checkpoint filename and `T`/`F` (update/don't update the file):

```sh
# Read density matrix from resid.dat (and rewrite it after)
./path-to-SSAdNDP/ssadndp.exe SSAdNDP.inp resid.dat

# Read from resid.dat without rewriting
./path-to-SSAdNDP/ssadndp.exe SSAdNDP.inp resid.dat F
```

### Input File — Spin Non-Polarized (ISPIN=1)

<div class="input-block" markdown="1">
<div class="input-block-label">SSAdNDP.inp</div>

```
GENERAL SEARCH INPUT FILE
Prepare files for visualization: T/F?
F
Number of different types of bonds in general search:
4
Number of centers, ON threshold, Distance threshold:
1 1.70 0
2 1.90 4.0
3 1.90 7.0
4 1.75 7.0
```

</div>

Line 3: whether to write visualization files (may be slow for large systems). Line 5: number of bond types. Lines 7+: for each bond type, specify *n*, ON threshold, and distance threshold (Å). In the example above, the program searches for lone pairs (ON &gt; 1.70), then 2c-2e, 3c-2e, and 4c-2e bonds in sequence.

### Input File — Spin Polarized (ISPIN=2)

<div class="input-block" markdown="1">
<div class="input-block-label">SSAdNDP.inp (spin-polarized)</div>

```
GENERAL SEARCH INPUT FILE
Prepare files for visualization: T/F?
F
Number of different types of bonds in general search (spin 1):
2
Number of centers, ON threshold, Distance threshold (spin 1):
1 0.8 0.
2 0.9 7.
Number of different types of bonds in general search (spin 2):
2
Number of centers, ON threshold, Distance threshold (spin 2):
1 0.8 0.
2 0.89 7.
```

</div>

For spin-polarized calculations, ONs are expected near 1.00 \|e\|. Different thresholds can be specified for each spin channel.

### Output and Checkpoint Files

An `SSAdNDP.out` file is written in the working directory. "Initial density" is the sum of diagonal elements of the current density matrix; if calculated from scratch, it equals the total number of electrons per unit cell. For each accepted bond, the ON and atom numbers (with cell indices relative to the central cell (0,0,0)) are written.

To save the initial density matrix to a checkpoint file for later trial runs:

```sh
# Write initial density matrix to checkpoint
./path-to-SSAdNDP/ssadndp.exe SSAdNDP.inp chkfilename

# Use it in a trial run without overwriting
./path-to-SSAdNDP/ssadndp.exe SSAdNDP.inp chkfilename F
```

<div class="note" markdown="1">
**Tip:** When resuming from a checkpoint, the "Initial density" at the top of the new `SSAdNDP.out` should equal the "Residual density" at the bottom of the previous run's output.
</div>
