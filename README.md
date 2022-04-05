# Solid State Adaptive Natural Density Partitioning (SSAdNDP)
*This project is originally written by Alexander Boldyrev @ Utah State University.
Forked from [🔗LINK](http://ion.chem.usu.edu/~boldyrev/ssadndp.php).*

*The source code has been updated by replacing several deprecated routine-calling scheme, and a modified Makefile for it to work with newest intel compiler and MKL version 2022.1.1*

Solid State Adaptive Natural Density Partitioning (SSAdNDP) [1] is an extension of the [AdNDP method](http://ion.chem.usu.edu/~boldyrev/adndp.php) [2] to periodic systems and as such was derived from [periodic implementation](http://schmidt.chem.wisc.edu/nbosoftware) [3] of the [Natural Bond Orbital (NBO)](http://nbo6.chem.wisc.edu/) analysis [4]. SSAdNDP allows the interpretation of chemical bonding in systems with translational symmetry in terms of classical lone pairs and two-center bonds, as well as multi-center delocalized bonding. Similar to AdNDP, the bonding pattern is expressed as a set of *n*-center - 2-electron (*n*c-2e) bonds. See the [AdNDP page](http://ion.chem.usu.edu/~boldyrev/adndp.php) and [AdNDP](http://pubs.rsc.org/en/Content/ArticleLanding/2008/CP/b804083d#!divAbstract) and [SSAdNDP](http://pubs.rsc.org/en/content/articlelanding/2013/cp/c3cp50350j#!divAbstract) publications [1,2] for details on this approach to interpretation of chemical bonding.

## Getting Started

These instructions will get you a copy of the project up and running.

This package comes in three part:
1. DFT interface and projection executable.
2. SSAdNDP calculation executable.
3. Visualization executable.

### Prerequisites

The code is written in `Fortran 90`, and requires the `BLAS95` and `LAPACK95` libraries as implemented in `MKL`, as well as libraries used for OpenMP parallelization.  

### Installing DFT interfaces
1. VASP (compatible with VASP v5.4.4)
  * copy `interfaces/VASP_interface/projection_output.F` into directory containing VASP source code `./scr`.
  * apply `interfaces/VASP_interface/proj_VASP_5.4.patch` by putting it in VASP source code's root directory and:
  ```
  patch -p0 < proj_VASP_5.4.patch
  ```
  * Compile VASP as usual.


2. CRYSTAL
  * No modification needed.

### Installing SSAdNDP software
1. Go to `SSAdNDP` folder, modify `Makefile` by changing the corresponding `MKLROOT`.
2. Use the following command to generate a `Makefile`:
```
./mkmf -t Makefile.template -p ssadndp.exe
```
3. Type `make` to get `ssadndp.exe`

4. To compile `projection.exe`, go to the `projection` folder and modify the file `Makefile` by changing the corresponding `MKLROOT`, then type:
```
Make
```

5. To compile `visual.exe`, go to the `visualization` folder and type:
```
ifort visual.f90 -o visual.exe
```



## Usage
The workflow of SSAdNDP is:
1. Calculate the system with DFT code and generate the wavefunctions (wavefunction.dat)
    1. Put `LNBO=.T.` into `INCAR`. Make sure `NPAR` is the default value and `ISYM=0` is set.
    2. Run VASP as usual and the code will generate `wavefunction.dat`.
2. Projecting the wavefunctions onto atomic orbitals (basis.inp) to get the projection coefficients (NBO.out)
    1. Construct a `basis.inp` using [basisexchange](https://www.basissetexchange.org) (note: choose wisely for this is very important and can have huge impact on later analysis). Remember to put `****` at the begining of the first species and make sure the ordering of species matches with the ordering of VASP's `POTCAR`.
    2. Run `projection.exe` by:
    ```
    projection.exe basis.inp wavefunction.dat NBO.out
    ```
3. Calculate the NBOs with predefined parameters (SSAdNDP*.inp)
    1. Construct a `SSAdNDP.inp` (or `SSAdNDP-UD.inp`) file for general search (or user-directed search). A set of exemplary inputs can be found in `manual`.
    2. Run SSAdNDP to generate `vis_*_bond.out` by (note `ulimit -s unlmited` maybe needed):
    ```
    ssadndp.exe SSAdNDP.inp
    ```
4. (optional) Plot the NBOs:
    1. run visualization program by:
    ```
    visual.exe vis_ud_bonds.out
    ```
    2. use vesta to visualize, first open the `lattice_vec.cube` and unselect `show section`, then `Edit` -> `Edit Data` -> `phase` -> `Import` to import `*_bond_*.cube`. Select the first entry of `*_bond_*.cube` and change the three numbers under `this layer` to `0.5 0.5 0.5`.

Please checkout `SSAdNDP_manual.pdf` under `manual`.

## References

The author request that anyone who downloads and utilizes the code cite:

[1] [Solid state adaptive natural density partitioning: a tool for deciphering multi-center bonding in periodic systems](http://pubs.rsc.org/en/content/articlelanding/2013/cp/c3cp50350j#!divAbstract) Timur R. Galeev, Benjamin D. Dunnington, J. R. Schmidt, and Alexander I. Boldyrev *Phys. Chem. Chem. Phys.*, 2013, **15**, 5022-5029

[2] [Developing paradigms of chemical bonding: adaptive natural density partitioning](http://pubs.rsc.org/en/Content/ArticleLanding/2008/CP/b804083d#!divAbstract) Dmitry Yu. Zubarev, and Alexander I. Boldyrev *Phys. Chem. Chem. Phys.*, 2008, **10**, 5207-5217

[3] [Generalization of Natural Bond Orbital Analysis to Periodic Systems: Applications to Solids and Surfaces via Plane-Wave Density Functional Theory](http://pubs.acs.org/doi/abs/10.1021/ct300002t) Benjamin D. Dunnington and J. R. Schmidt *J. Chem. Theory Comput.*, 2012, **8**, 1902-1911

[4] See [NBO references and bibliography](http://www.chem.wisc.edu/~nbo5/biblio.htm)

## License

  This project is licensed under the GNU License - see the `LICENSE.md` for details
