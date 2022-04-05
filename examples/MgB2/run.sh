
mpirun -n 10 /home/cx219/source/VASP/2022-03-31/vasp.5.4.4_NBO/bin/vasp_std

./bin/projection basis.inp wavefunction.dat NBO.out
./bin/ssadndp.exe SSAdNDP-UD.inp
./bin/visual.exe vis_ud_bonds.out
