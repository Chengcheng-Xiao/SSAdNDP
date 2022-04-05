base=projection.exe
compiler=ifort
MKLROOT=/share/compilers/intel/oneapi/mkl/2022.0.1/

projection_main:  projection_main.f90
	${compiler} rd_wavefunction.f90  projection_shared.f90 rd_basis.f90 bloch_overlap.f90 PAW.f90 projection_main.f90  -o ${base} -I${MKLROOT}/include/intel64/lp64/ -L${MKLROOT}/lib/intel64 -lmkl_lapack95_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread -lmkl_blas95_lp64 -qopenmp



