Example 1 - MgB2

This example is based on [Phys. Chem. Chem. Phys., 2013,15, 5022-5029](https://pubs.rsc.org/en/content/articlelanding/2013/cp/c3cp50350j)

Steps to reproduce the results:
1. Run VASP.
2. Run productions:
```
projection.exe basis.inp wavefunction.dat NBO.out
```
3. Run SSAdNDP (general search):
```
ssadndp.exe SSAdNDP.inp
```
4. Run SSAdNDP (User-directed):
```
ssadndp.exe SSAdNDP-UD.inp
```
5. Get visualization:
```
visual.exe vis_ud_bonds.out
```
