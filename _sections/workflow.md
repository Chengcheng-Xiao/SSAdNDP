---
anchor: workflow
title: Workflow
order: 3
---

## Workflow

<ol class="steps">
<li>
<div class="step-body" markdown="1">
**DFT calculation — generate wavefunctions**

In your VASP `INCAR`, set `LNBO=.T.`, ensure `NPAR` is at its default, and set `ISYM=0`. Run VASP; it will produce `wavefunction.dat`.
</div>
</li>
<li>
<div class="step-body" markdown="1">
**Project wavefunctions onto atomic orbitals**

Construct `basis.inp` using [Basis Set Exchange](https://www.basissetexchange.org). Put `****` at the start of the first species; match species ordering to VASP's `POTCAR`.

```
projection.exe basis.inp wavefunction.dat NBO.out
```

</div>
</li>
<li>
<div class="step-body" markdown="1">
**Run SSAdNDP analysis**

Construct `SSAdNDP.inp` (general search) or `SSAdNDP-UD.inp` (user-directed). See sections below for input formats.

```
ssadndp.exe SSAdNDP.inp
```

<div class="note" markdown="1">
**Note:** `ulimit -s unlimited` may be needed when running `ssadndp.exe` on large systems.
</div>
</div>
</li>
<li>
<div class="step-body" markdown="1">
**Visualize (optional)**

```
visual.exe vis_ud_bonds.out
```

Then open in VESTA — see the [Visualization](#visualization) section.
</div>
</li>
</ol>
