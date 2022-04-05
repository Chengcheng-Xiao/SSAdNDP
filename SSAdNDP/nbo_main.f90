!!! TRG: this is a shortened version of BDD's nbo_main.f90 code
!!! It uses the following BDD's procedures: 
!!! -read_input_file, 
!!! -calc_nelec, 
!!! -do_pre_nao, 
!!! -do_nao 
!!! to obtain the density matrix in NAO basis for the SSAdNDP analysis.
!!! Then it executes the SSAdNDP analysis instead of running the NBO analysis.
!!!  few modifications made to this part of code are marked with 
!!! '!!! TRG: ...'



!This is an NBO analysis code appropriate for periodic systems

PROGRAM nbo_main
  USE periodic_matutil
  USE nbo_shared
  USE matutil
  USE pre_nao
  USE nao
  USE nbo
  !!! TRG:
  USE ssadndp_module      



  IMPLICIT NONE

  !Command line arguments
  CHARACTER*128 :: buffer, inputfile !!! TRG: added inputfile for SSAdNDP.inp or SSAdNDP-UD.inp

  !Define total number of atoms, basis functions, and (relevant) real-space matrices (i.e. overlap and density)
  !See PRB, 61, 16440, "Linear-scaling DFT with Gaussian orbitals...", Kudin and Scuseria for details of the
  !real-space implementation of periodic SCF and definitions of matrices S^{og}_{mu nu} and P^{og}_{mu nu}
  INTEGER :: natom,nbasis,ng,nk,nspins
  INTEGER,ALLOCATABLE ::  nnbo(:)
  !The mapping of the basis functions to atoms (i.e. first basis function for each atom)
  !as well as the angular momentum quantum number (l) for each basis function
  INTEGER,DIMENSION(:),ALLOCATABLE,TARGET :: ibasismap,ishellmap,ilmap,immap,iatnum
  REAL*8,DIMENSION(:),ALLOCATABLE,TARGET  :: iatval
  !The mapping of the ng matrices into spatial cells
  INTEGER,DIMENSION(:,:),ALLOCATABLE,TARGET ::  indexg
  !The k-points to used, as coefficients of reciprocal lattice vectors
  REAL*8,DIMENSION(:,:),ALLOCATABLE,TARGET  ::  kpt
  !The weight of each k-point as used in the plane wave calculation.  From space group of the Brillouin zone
  REAL*8,DIMENSION(:),ALLOCATABLE,TARGET    ::  kpt_wt
  !A listing of the symbol of each atom
  CHARACTER*2,DIMENSION(:),ALLOCATABLE,TARGET :: symbols
  !Will hold the calculated real-space density and overlap matrices, rho^{0g}_{mu nu}
  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE,TARGET :: rho0, fock0
  REAL*8,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: s0
  !The initial k-space density and overlap matrices
  COMPLEX*16,DIMENSION(:,:,:,:),ALLOCATABLE,TARGET :: rhok, fockk
  COMPLEX*16,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: sk
  REAL*8,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: transform  !Will hold the real space transformation matrix to go from NHO's back to AO basis

  REAL*8,ALLOCATABLE      ::   output_coeff(:,:,:)  !Hold the coefficients of each orbital in the AO basis, can then be used for visualization
  COMPLEX*16,ALLOCATABLE  ::   bloch_coeff(:,:,:,:)
  COMPLEX*16,ALLOCATABLE  ::   fock_nbo(:,:,:,:), rho_nbo(:,:,:,:)
  REAL*8,ALLOCATABLE      ::   rho_dummy(:,:,:,:), real_fock_nbo(:,:,:,:)

  REAL*8,ALLOCATABLE    ::  output_occ(:)
  REAL*8                ::  energy_diff, occ_diff,perturb,occ_factor

  CHARACTER(32),ALLOCATABLE  ::   nbo_lbl(:)

  COMPLEX*16,ALLOCATABLE    ::   energy_test(:,:)
  COMPLEX*16                ::   energy_sum, occ_sum  !For testing bloch-space matrices in NBO basis

  REAL*8     :: ti,tf

  !Checkpoint file information
  LOGICAL :: checkpoint_exists, write_checkpoint

  !A structure which encapsulates all that information
  TYPE(nbo_input) :: inp

  !Temporary variable
  INTEGER :: ig,inbo,ik,i,j,k,ispin
  COMPLEX*16  ::  arg

  INTEGER,ALLOCATABLE  :: nbond(:),nlp(:),nryd(:)

  REAL*8,PARAMETER           ::  pi=4.d0*ATAN(1.d0)
  COMPLEX*16,PARAMETER       ::  sqrt_minus_one=(0.d0, 1.d0)
  REAL*8,PARAMETER           ::  aukcal=627.509469d0


!!!TRG: allowing third command-line argument to specify whether to update the
!resid matrix file:
  character*1::arg3
  integer:: istat, istat1
   

  CALL CPU_TIME(ti)



!!! TRG: making the nbo input file default to NBO.out:
  open (10, file='NBO.out', status='old')
!!!  !Do IO; first read the input file
!!!  CALL GETARG(1,buffer)
!!!  OPEN(10, FILE=buffer, STATUS='old')
  CALL read_input_file
  CLOSE(10)

!!! TRG: read the SSAdNDP.inp or SSAdNDP-UD.inp file name as the first command-line argument; those will be used for general or user-directed search
  call getarg (1, inputfile)
  if ((inputfile.ne.'SSAdNDP.inp').and.(inputfile.ne.'SSAdNDP-UD.inp')) then
    write (*,*) 'ERROR: smth wrong with SSAdNDP input-file name?'
    write (*,*) 'Stopping'
    STOP
  end if




  !See if a checkpoint file was specified

!!!: TRG slightly modifying the file handling and
!allowing third command-line argument to specify whether to update the
!resid matrix file:
   checkpoint_exists=.false.

   IF (IARGC().GT.1) THEN
     CALL GETARG(2,buffer)
!!!TRG     OPEN(10, FILE=buffer, STATUS='old', FORM='unformatted', ERR=20)

     write (*,*)
     open(10, file=buffer, status='old', form='unformatted', iostat=istat)
     if (istat==0) then
       write(*,*) 'Reading (residual) density matrix in NAO basis from checkpoint file: ', buffer
!!!TRG       READ(10) inp%rho0
!!!TRG       READ(10) inp%transform
!!!TRG       READ(10) inp%fock0
       read(10, iostat=istat1) inp%rho0
       if (istat1.ne.0) then
         write(*,*) 'Error reading from', buffer
         write (*,*) 'Stopping'
         stop
       end if
       read(10, iostat=istat1) inp%transform
       if (istat1.ne.0) then
         write(*,*) 'Error reading from', buffer
         write (*,*) 'Stopping'
         stop
       end if
       read(10, iostat=istat1) inp%fock0 
       if (istat1.ne.0) then
         write(*,*) 'Error reading from', buffer
         write (*,*) 'Stopping'
         stop
       end if

       inp%s0=periodic_matiden(inp%nbasis)
       checkpoint_exists=.true.
     else
       write (*,'(1X, A33, A30)') 'Looks like checkpoint file named ', buffer
       write (*,*) "doesn't exist"
       write (*,*) 'Density matrix in NAO basis will be calculated from scratch'
       write (*,*)
       checkpoint_exists=.false.
     end if
       
     close(10)

     write_checkpoint=.true.

     if (iargc()>2) then
       call getarg (3, arg3)
       read (arg3,*) write_checkpoint
     end if

 
     if (write_checkpoint) then
       write (*,'(1X, A57, A30)') 'The (updated) residual density matrix will be written to ', buffer
     else 
       write (*,*) 'The checkpoint file will not be updated'
     end if
     write (*,*)
   ENDIF

  !Check number of electrons and do transform to an orthogonal basis
  CALL calc_nelec(inp,checkpoint_exists)

  IF (.NOT.checkpoint_exists) THEN
     !Convert to a symmetrically orthoganlized basis, weighted by the occupations of the pre-NAOs.
     !This is the NAO basis.
     CALL do_pre_nao(inp)
     CALL do_nao(inp)           !!! TRG: moved this line here for SSAdNDP
  ENDIF

  !!!TRG:
  call SSAdNDP(inp, nspins, inputfile) 

  IF (write_checkpoint) THEN
     write(*,*)
     write(*,*) 'Writing (residual) density matrix in NAO basis to checkpoint file: ', buffer

     OPEN(10, FILE=buffer, STATUS='replace', FORM='unformatted')

     WRITE(10) inp%rho0
     WRITE(10) inp%transform
     WRITE(10) inp%fock0
  else 
     write(*,*)
     write (*,*) 'The current NAO density matrix will not be saved'
     write (*,*)
  ENDIF
  CLOSE(10)

  !!!TRG:
  stop   
                       




     



  !Now get down to buisiness; first output the NAOs (this is easy, since we are already in the
  !NAO basis due to the prior transformations
  CALL do_nao(inp)
  
  !***This has been removed so that all one center orbitals are treated as lone pairs.
  !***The lone pairs were not actually projected, but acutally depleted. Less rigorous
  !***It was actually not even depletion (which would be equivalent in this case) just zeroing diagonals
  !Project off the core electrons (with occupancy ~2) from the density matrix
  !CALL do_core_projection(inp)
  

  ALLOCATE(rho_dummy(nbasis,nbasis,ng,nspins))
  rho_dummy = inp%rho0

  !Now do the final NBO analysis
  !ALLOCATE(bond_out(inp%nspins),lp_out(inp%nspins),ryd_out(inp%nspins))
  ALLOCATE(output(nspins),nnbo(nspins))
  DO ispin=1,nspins
     CALL do_nbo(inp,ispin,nnbo(ispin))
  ENDDO
  CALL CPU_TIME(tf)
  WRITE(6,*)'total time for NBO analysis',SNGL(tf-ti)


  !!!!!!!!!!!!!!!!!
  !This is STOP #1!
  !!!!!!!!!!!!!!!!!
  STOP


  !This section of the code is for visualization purposes.

  !Convert the nbo orbitals with coeff in the NAO basis into coeff in the AO basis for visualization.
  DO ispin=1,nspins
     DO inbo=1,nnbo(ispin)
        output(ispin)%coeff(:,inbo,:) = periodic_matvecmul(transform,output(ispin)%coeff(:,inbo,:))
     ENDDO
  ENDDO

  !Now actually write out a file to be read in for use in visualization
  OPEN(95,file='nbo_vis.out')
  WRITE(95,*)"Output of lattice vector and ao coeffs from JRS periodic NBO code"
  WRITE(95,*)

  !System information
  WRITE(95,*)nbasis, '! number of basis functions per unti cell'
  WRITE(95,*)nnbo(1), '! number of popssible lonepairs and NBOs, set of coefficients'
  WRITE(95,*)ng, '! number of l_vectors, unit cell pairs'
  WRITE(95,*)nspins,'! number of unique spins for which NBOs have been obtained'
  WRITE(95,*)

  !Real space unit cells
  DO ig=1,ng
     WRITE(95,*)indexg(:,ig)
  ENDDO
  WRITE(95,*)

  !Coefficients of each NBO in the AO basis
  DO ispin=1,nspins
     DO ig=1,ng
        DO inbo=1,nnbo(1)
           WRITE(95,*)output(ispin)%coeff(:,inbo,ig)
        ENDDO
        WRITE(95,*)
     ENDDO
  ENDDO

  CLOSE(95)


  !!!!!!!!!!!!!!!!!
  !This is STOP #2!
  !!!!!!!!!!!!!!!!!
  STOP


  !This section of the code is for energy analysis in the NBO basis, which is currently not numerically meaningful.
  !The only useful thing is that the sum of the trace of [rho*fock] summed over k-points should still match the original VASP value.

  !Note as this is currently formatted, it is not compatible with the visualization output above.
  !The above process destructively places, the NBO coeffs in the AO basis in 'output'.
  !The below process assumes the coeffs are in the NAO basis.

  !Check to make the sure a full number of NBO's was found. 
  !Basically did a bond get rejected after orthogonalization.
  !A square matrix is necessary for the unitary transforms.
  DO ispin=1,nspins
     IF( nnbo(ispin) /= nbasis )THEN
        WRITE(6,*)'Incorrect number of NBOs for spin',ispin
     ENDIF
  ENDDO

  ALLOCATE(rho_nbo(nbasis,nbasis,nk,nspins),fock_nbo(nbasis,nbasis,nk,nspins),real_fock_nbo(nbasis,nbasis,ng,nspins))
  ALLOCATE(bloch_coeff(nbasis,nbasis,nk,nspins))


  inp%fockk = 0.d0
  inp%rhok = 0.d0
  bloch_coeff = 0.d0

  IF( real_init )STOP 'Real space intialized calculations are not setup for past bond search'

  !Go to Bloch space (unitary transform will then be straight forward)
  DO ispin=1,nspins
     CALL real_to_bloch(inp%fockk(:,:,:,ispin),inp%fock0(:,:,:,ispin),inp%kpt,inp%indexg)
     CALL real_to_bloch(inp%rhok(:,:,:,ispin),rho_dummy(:,:,:,ispin),inp%kpt,inp%indexg)
     CALL real_to_bloch(bloch_coeff(:,:,:,ispin),output(ispin)%coeff,inp%kpt,inp%indexg)
  ENDDO

  ALLOCATE(energy_test(nbasis,nbasis))
  energy_sum = 0.d0
  occ_sum = 0.d0

  !Perform unitary transform and multiplications to test sum of orbital energies
  DO ik=1,nk
     DO ispin=1,nspins
        fock_nbo(:,:,ik,ispin)=matunitary_trans(inp%fockk(:,:,ik,ispin),TRANSPOSE(bloch_coeff(:,:,ik,ispin)))
        rho_nbo(:,:,ik,ispin)=matunitary_trans(inp%rhok(:,:,ik,ispin),TRANSPOSE(bloch_coeff(:,:,ik,ispin)))

        energy_test = MATMUL(fock_nbo(:,:,ik,ispin),rho_nbo(:,:,ik,ispin))
        DO j=1,nnbo(ispin)
           energy_sum = energy_sum + energy_test(j,j)*inp%kpt_wt(ik)
           occ_sum = occ_Sum + rho_nbo(j,j,ik,ispin)*inp%kpt_wt(ik)
        ENDDO
     ENDDO
  ENDDO
  WRITE(6,*)'energy test in nbo basis set    ',energy_sum  
  WRITE(6,*)'occupancy test in nbo basis set ',occ_sum
  WRITE(6,*)

  !Then take matrixces in NBO basis back to real space
  DO ispin=1,nspins
     CALL bloch_to_real(inp,fock_nbo(:,:,:,ispin),real_fock_nbo(:,:,:,ispin),inp%kpt,inp%indexg)
     CALL bloch_to_real(inp,rho_nbo(:,:,:,ispin),rho_dummy(:,:,:,ispin),inp%kpt,inp%indexg)
  ENDDO

  !Dump out some information on the NBO's (diagonal amtrix elements)
  DO ispin=1,nspins
  IF( nspins .GT. 1 )WRITE(6,*)"NBOs for spin",ispin
  WRITE(6,'(A,4I5)')'fock matrix of nbos in central unit cell'
  DO inbo=1,nnbo(ispin)
     WRITE(6,*)inbo,real_fock_nbo(inbo,inbo,1,ispin)
  ENDDO
  WRITE(6,*)
  
  WRITE(6,*)'density matrix of nbos in central unit cell'
  DO inbo=1,nnbo(ispin)
     WRITE(6,*)inbo,rho_dummy(inbo,inbo,1,ispin)
  ENDDO
  WRITE(6,*)
  ENDDO


  !!!!!!!!!!!!!!!!!
  !This is STOP #3!
  !!!!!!!!!!!!!!!!!
  STOP


  !This is second order perturbation analysis.
  !This suffers from very high numerical inaccuracies based on the projection process.
  !From testing, the results seem to be what you would expect
  !i.e. back bonding from metal surfaces into pi* orbitals being the largest factor
  !But I don't trust any of it.  
  !This relies on the real space Fock matrices obtained in the above step

  DO ispin=1,nspins

  IF( nspins .GT. 1 )WRITE(6,*)'perturbation analysis for spin',ispin

  !Looking at perturbations driven donations from j -> i nbo 
  DO j=1,nnbo(ispin)

     !Make sure the orbital is capable of donation, based on a high occupancy
     IF( output(ispin)%occ(j) < 0.7d0*DBLE(3-nspins) )GOTO 9 

     DO i=1,nnbo(ispin)

        !Check that the acceptor orbital is unique and has a low enough occupancy to accept
        IF( i == j )GOTO 40
        IF( output(ispin)%occ(i) > 0.7d0*DBLE(3-nspins) )GOTO 40

        energy_diff = real_fock_nbo(i,i,1,ispin) - real_fock_nbo(j,j,1,ispin)
        occ_diff = output(ispin)%occ(j) - output(ispin)%occ(i)

        !WRITE(6,'(A,I4,A15,I4,A15)')'looking at nbos',j,output(ispin)%label(j),i,output(ispin)%label(i)
        !WRITE(6,*)'donation is possible'
        !WRITE(6,*)
        !Since donation is possible, loop over all possible unit cells to see if a big enough delocalization exists
        DO ig=1,ng

           !The occupancy used is NOT simply the occupancy of the donor orbital
           !The perturbative energy lowering comes from a mixing of the occupancy of the two orbitals, thus their occupancies should be added
           !However if it will result in an occupancy over 2, some electrons will go to the higher energy split and cancel out the lowering effect
           !Thus occ_factor represents how many electrons will be serving to lower the energy in the resulting orbital mixing
           occ_factor = output(ispin)%occ(j) + output(ispin)%occ(i)
           IF( occ_factor > DBLE(3-nspins) )occ_factor = 2.d0*DBLE(3-nspins) - occ_factor
           perturb = occ_factor * real_fock_nbo(j,i,ig,ispin)**2 / energy_diff
           perturb = perturb * aukcal

           IF( perturb > 0.25d0*DBLE(3-nspins) )THEN
               WRITE(6,'(A,I5,A15,A,I5,A15,I4)')'Relevant perturbation found ',j,output(ispin)%label(j),'to ',i,output(ispin)%label(i),ig
               WRITE(6,*)'Energy diff and coupling',SNGL(energy_diff),SNGL(real_fock_nbo(j,i,ig,ispin))
               WRITE(6,*)'pertubation energy',SNGL(perturb),'kcal/mol'
               WRITE(6,*)
           ENDIF

        ENDDO

70      CONTINUE

        !WRITE(6,*)

40   ENDDO
9 ENDDO


  ENDDO


  STOP


CONTAINS

  !
  !Reads the input file and allocates necessary memory
  !All matrices read in are in reciprocal space, NOT real space
  !
  SUBROUTINE read_input_file
    USE matutil
    IMPLICIT NONE
    INTEGER :: ibasis,im,ik,ispin
    !LOGICAL :: real_init
    !INTEGER :: nkx,nky,nkz
    REAL*8  :: num_elec

    CHARACTER(64)  :: NBO_mat_fn

      real_init = .TRUE.

      READ(10,*) !read comment line
      READ(10,*) natom
      READ(10,*) nbasis
      READ(10,*) nspins
      READ(10,*) ng

      ALLOCATE(ibasismap(natom+1))
      READ(10,*,ERR=20)ibasismap
      GOTO 30

!20    WRITE(6,*)'The input file has an entry for nk, so the matrices must be in bloch space'
20    CONTINUE
      real_init = .FALSE.
      nk = ibasismap(1)


30    CONTINUE
      !WRITE(6,*)'real start boolean   ',real_init


      IF( real_init )STOP 'NBO code is no longer compatible with Gaussian output'

      ALLOCATE(ishellmap(nbasis))
      ALLOCATE(ilmap(nbasis))
      ALLOCATE(immap(nbasis))
      ALLOCATE(symbols(natom))
      ALLOCATE(iatnum(natom))
      ALLOCATE(iatval(natom))
      ALLOCATE(indexg(3,ng))

      IF( .NOT. real_init )THEN
         ALLOCATE(kpt(3,nk))
         ALLOCATE(kpt_wt(nk))
         READ(10,*) ibasismap
      ENDIF
      READ(10,*) ishellmap
      READ(10,*) ilmap
      READ(10,*) symbols
      READ(10,*) iatnum
      IF( real_init )THEN
         iatval = iatnum
      ELSE
         READ(10,*) iatval
      ENDIF

      !WRITE(6,*)'iatval',iatval

      !WRITE(6,*)'made it to reading g-vectors'

      !Read in information about real space g-vectors.
      !These will only be for nearest neighbor cells where bonds we will be searched for
      DO ig =1,ng
         READ(10,*) indexg(1,ig), indexg(2,ig), indexg(3,ig)
      ENDDO
      CALL periodic_matinit(ng,indexg)

      IF( real_init )THEN
         CALL get_nk(nkx,nky,nkz)
         nk=nkx*nky*nkz
      ENDIF


      ALLOCATE(rho0(nbasis,nbasis,ng,nspins))
      ALLOCATE(s0(nbasis,nbasis,ng))
      ALLOCATE(fock0(nbasis,nbasis,ng,nspins))
      ALLOCATE(transform(nbasis,nbasis,ng))
      
      ALLOCATE(rhok(nbasis,nbasis,nk,nspins))
      ALLOCATE(sk(nbasis,nbasis,nk))
      ALLOCATE(fockk(nbasis,nbasis,nk,nspins))

      !WRITE(6,*)'allocated everything'

      IF( real_init )THEN

          STOP 'Need to resetup matrix input for real space initialization systems'


!         CALL read_real_triangular(10,nbasis,s0(:,:,1),s0(:,:,1))
!         DO ig=2,ng,2
!            CALL read_real_triangular(10,nbasis,s0(:,:,ig),s0(:,:,ig+1))
!            CALL read_real_triangular(10,nbasis,s0(:,:,ig+1),s0(:,:,ig))
!         ENDDO
!         CALL read_real_triangular(10,nbasis,rho0(:,:,1),rho0(:,:,1))
!         DO ig=2,ng,2
!            CALL read_real_triangular(10,nbasis,rho0(:,:,ig),rho0(:,:,ig+1))
!            CALL read_real_triangular(10,nbasis,rho0(:,:,ig+1),rho0(:,:,ig))
!         ENDDO
!
!         OPEN(12,file='NBODATA.54')
!         DO j=1,3
!            READ(12,*)
!         ENDDO
!         fock0=0.d0
!         CALL read_real_triangular(12,nbasis,fock0(:,:,1),fock0(:,:,1))
!         !WRITE(6,*)'test of read in fock matrix'
!         !DO ibasis=1,nbasis
!         !   WRITE(6,'(27F10.5)')fock0(ibasis,:,1)
!         !   WRITE(6,*)
!         !ENDDO
!         
!         CLOSE(12)
!
!         !STOP 'Fock matrix input has not been implemented for real space start'
!
!         DO ig=1,ng
!            num_elec=num_elec + mattrace(MATMUL(rho0(:,:,ig),TRANSPOSE(s0(:,:,ig))))
!         ENDDO
!         WRITE(6,*)'number of initial electrons in real space',2.d0*num_elec
!         num_elec = 0.d0
!         DO ig=1,ng
!            num_elec=num_elec + mattrace(MATMUL(rho0(:,:,ig),TRANSPOSE(fock0(:,:,ig))))
!         ENDDO
!         WRITE(6,*)'initial real space energy sum',2.d0*num_elec
!
!         sk=0.d0
!         rhok=0.d0
!         sk=periodic_matbloch(s0,nkx,nky,nkz)
!         rhok=periodic_matbloch(rho0,nkx,nky,nkz)
!         fockk=periodic_matbloch(fock0,nkx,nky,nkz)

      ELSE

         !WRITE(6,*)'reading in k-pts'

         !Read in k-vectors indexing the overlap and density matrices.
         !These are in the same order as the matrices 
         DO ik=1,nk
            READ(10,*)kpt(:,ik),kpt_wt(ik)
            !WRITE(6,*)ik,kpt_wt(ik)
         ENDDO

         READ(10,*)NBO_mat_fn
         CLOSE(10)

         !WRITE(6,*)'Input matrices will be read in from',NBO_mat_fn

         OPEN(66,file=NBO_mat_fn,FORM='UNFORMATTED')

         !read in the overlap and density matrices s^{k}_{mu nu}
         !At each k-point each of these matrices are hermitian
         !Therefore only half of each matrix is written out in the file, and symmetry is used to file in the remainder
         DO ik=1,nk
            CALL read_bloch_triangular(66,nbasis,sk(:,:,ik))
         ENDDO

         DO ispin=1,nspins
            DO ik=1,nk
               CALL read_bloch_triangular(66,nbasis,rhok(:,:,ik,ispin))
            ENDDO
         ENDDO

         DO ispin=1,nspins
            DO ik=1,nk
               CALL read_bloch_triangular(66,nbasis,fockk(:,:,ik,ispin))
            ENDDO
         ENDDO

      ENDIF

      !correct for the fact we (maybe) just read the alpha density
      rhok=rhok * DBLE(3-nspins)
   
      !Finally, initalize immap to contain the m quantum number for each basis function
      !Assumes first basis function is s-type
      immap=0
      DO ibasis=2,nbasis
         IF (ilmap(ibasis).NE.ilmap(ibasis-1).OR.ishellmap(ibasis).NE.ishellmap(ibasis-1)) THEN
            DO im=0,2*ilmap(ibasis)
               immap(ibasis+im)=im
            ENDDO
         ENDIF
      ENDDO

      !Now store these values in our "system" structure to keep them all together
      inp%natom=natom
      inp%nbasis=nbasis
      inp%ng=ng
      inp%nk=nk
      inp%nspins=nspins
      inp%ibasismap=>ibasismap
      inp%ishellmap=>ishellmap
      inp%ilmap=>ilmap
      inp%immap=>immap
      inp%indexg=>indexg
      IF( .NOT. real_init )THEN
          inp%kpt=>kpt
          inp%kpt_wt=>kpt_wt
      ENDIF
      inp%symbols=>symbols
      inp%iatnum=>iatnum
      inp%iatval=>iatval
      inp%rho0=>rho0
      inp%s0=>s0
      inp%fock0=>fock0
      inp%transform=>transform
      inp%sk=>sk
      inp%rhok=>rhok
      inp%fockk=>fockk

    END SUBROUTINE read_input_file

END PROGRAM nbo_main


