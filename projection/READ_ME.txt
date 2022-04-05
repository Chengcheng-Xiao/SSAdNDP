***********************************************************
***********************************************************

Projection Code Instructions and guide lines
By BDD 2/1/12

11/12/12 - Update output file description

A description of the projection can be found in:
Dunnington and Schmidt, JCTC; 8, 1902-1911 (2012)
***********************************************************
***********************************************************

This projection code is capable of representing any plane-wave bands output by VASP into 
any user specified Gaussian type orbital basis set.
The code is currently set up to then output all information necessary for periodic NBO 
analysis using the code of JRS and BDD.

VASP must be modified to output the necessary information that can then be
read by the projection code.  The included patch should be applied to fileio.F
to accomplish this.  VASP will now output and unformatted file called
wavefunction.dat containing all necessary information.  Since this file is
large output of it is optional.  To get the code to output the file, the
enviromental variable 'NBO_VASP' must be set to 'yes', via:
   export NBO_VASP=yes

The makefile contained with this file will produce an executable called 'projection'.  
This executable then runs with three (optional) inputs. 
They are listed below in order they will be read in by the program, along with the default 
file name the program will look for.

1. Basis set file. 'basis.inp' 
   This should contain all information on the atomic orbital basis set to be used, in the 
   format of Gaussian94. This means any header can be used, but each atom type must begin and end
   with a line of '****'.
   Basis set orbital types can include everything up to f-type (including 'sp'-labeling). 
   Within each orbital, it is assumed that the individual Gaussians are listed in
   descending order. This is used for screening purposes in the following subroutines:
      -real_overlap
      -nu_g_overlap
      -stype_PAW_offsite
      -ptype_PAW_offiste

2. VASP output file. 'wavefunction.dat'
   This file is output from our customized version of VASP as wavefunction.dat and 
   contains all necessary information about the plane-wave output for use in the projection.
   In running VASP there are a couple of limitations.

      i. ISYM=0, This parameter must be specified in the INCAR file.  It limits the symmetry 
         used in sampling the Brillouin zone to only inversion symmetry.  This is currently 
         all the projection code is capable of handling.

     ii. NO ULTRASOFT PSEUDOPOTENTIALS.  There is no functional form associated with this
         pseudopotential and it is therefore impossible to project a correction for the 
         valence electrons' interaction with the core electrons.
         This means the valence bands will each reprensent a (different) incorrect number of 
         electrons.
         Norm conserving and PAW type pseudopotentials are compatible with the code.
         However, with PAW pseudopotentials approximations are introduced that make it possible
         for the norm of projected bands to extend beyond 1.

    iii. Gamma-point containing k-point mesh. To take advantage of the inversion symmetry of the 
         Brillouin zone, I assume that the first k-point is the gamma point.

   Outside of this any type of calculation that can be performed with VASP is supported, 
   including spin polarized calculations.

3. NBO output file. 'NBO.out'
   This file contains information, in a readable format, necessary for the periodic NBO 
   code of JRS and BDD.
   The header of the file contains information on the basis set as well as atoms of the system.
   It also contains the k-points utilized in the VASP calculation (where the projection has 
   been performed) as well as the indices of the eral space unit cells that will be used in the
   NBO computation. This set of unit cells is only those which are nearset neighbors to the 
   central unit cell.
   The actual matrices are stored in 'NBO_mat.out'. For each k-point there is an overlap, 
   density and fock matrix (two density and fock for spin polarized calculations). 
   This file is unformatted for memory considerations as well as efficient input/output.
   The file name of the matrix storage is placed at the end of 'NBO.out' and read from there
   by the NBO code, so it will have to be changed if the file name is changed.


Output files:

The projection code outputs a couple additional files.

1.  'spillover.out' Contains quantitative information on the quality of the projection. 
     SPILLOVER is defined as in Eq 10 of the periodic NBO paper. For norm-conserving 
     pseudopotentials, this value is rigorously bound by 0 and 1. 
     However, for PAW-type pseudopotentials, approximations are made in calculating atomic 
     orbital-augmenter overlap and this constraint is lifted. As a result the spillover sum 
     includes both positive and negative terms. To account for this we have defined the SPREAD
     parameter, which is defined the same as spillover, except using the absolute value of all 
     terms within the summation.
     Both spillover and spread are calculated across all bands, as well as only those that are 
     occupied, and therefore contribute to the density matrix.  
     Additionally, NO norm should be above 1.01 and a tally of occupied bands who breach this
     limit is included in this file.

     Additionally, atomic weighted spread and spillover are provided based on occupied bands.
     These are calculated similarly to Eq 10, except there is now a weighting (and normalization 
     factor) calculated as the sum of the squares of the coefficients in the projected band of all 
     basis functions centered on a particular atom.

     Finally is a listing of the occupied band that the AO basis does the worst job in 
     representing, gauged by the spread of that band. These give an upward bound on the error 
     incurred using a projection


2.  'band_spillover.out' This contains the norm of each projected band, for each spin at each 
     k-point.
     This just gives a more detailed picture of the information summarized in the spillover.out 
     file.

The program also outputs information to the screen.
This mainly consists of system information from the VASP wavefunction file.
Additionally after the projection the density and fock matrices (in the projected matrices) are 
checked by calculating the associated observable, number of electrons and sum of orbital 
eigenvalues respectively. These are also calculated from the VASP input at the beginning of the 
calculation.






