PROGRAM visual
  IMPLICIT NONE

  TYPE AO_function
     REAL*8, ALLOCATABLE, DIMENSION(:) :: norm, alpha, coeff
     INTEGER                :: num_gauss
     !INTEGER, DIMENSION(3)  :: lmn
     INTEGER                :: atom
     REAL*8, DIMENSION(3)   :: pos
     INTEGER                :: level  !Keeps track of what basis functions share the same sets of exponents, on the same atom.
     INTEGER                :: l,m     !l- and m-quantum numbers of basis function.  m actually runs from 0,2*l+1 and is onyl and index.

     INTEGER                :: ncart
     INTEGER,ALLOCATABLE    :: cart_mat(:,:)
     REAL*8,ALLOCATABLE     :: cart_coeff(:)

  END TYPE AO_function

  TYPE(AO_function),ALLOCATABLE  ::  AO_basis(:)


  REAL*8,DIMENSION(3,3)  ::  a

  REAL*8,ALLOCATABLE     ::  ao_coeff(:,:,:,:)
  INTEGER, ALLOCATABLE   ::  index_l(:,:)

  INTEGER                ::  s_dim, ao_nbasis
  INTEGER                ::  nnbo, num_l
  INTEGER                ::  nspins

  INTEGER                ::  n_atom
  REAL*8, ALLOCATABLE    ::  atom_pos(:,:)
  INTEGER, ALLOCATABLE   ::  atomic_number(:)
  REAL*8, DIMENSION(3)   ::  image_pos
  REAL*8, DIMENSION(3,2) ::  image_screen
  INTEGER                ::  image_count, image_tot
  LOGICAL, ALLOCATABLE   ::  image_disp(:)

  CHARACTER(64)          ::  comment

  INTEGER                :: il,inbo,nu,ispin


  !The following variables are all used in calculating (numerically) and writing out the density of each band   
  !This can then be plotted to look at the shape of the resulting bands in space.
  !This can also be used to test for the norm of each band by numerical integration of the density
  REAL*8,ALLOCATABLE      ::  norm_test(:,:)
  REAL*8                  ::  gridvol
  CHARACTER(30),ALLOCATABLE ::  file_names(:,:)
  INTEGER                 ::  file_num
  INTEGER, DIMENSION(3)   ::  mesh  !mesh size to use
  REAL*8, DIMENSION(3,3)  ::  box   !box size to grid
  REAL*8, DIMENSION(3)    ::  r,origin!     !radius used in griding and origin for the cube file
  REAL*8                  ::  density
  REAL*8                  ::  gaussian, gauss_coeff
  REAL*8                  ::  rsqr
  REAL*8                  ::  cartesian
  REAL*8                  ::  wave_funct
  REAL*8,DIMENSION(3)     ::  r_pos
  REAL*8                  ::  screen
  INTEGER                 :: ix,iy,iz, i,j,k

!!!TRG Adding command line argument variables for filename and box multiplier
  character(28) :: inputfile
  character(4) :: arg2, arg3, arg4, arg5
  logical :: ifperiodic




  !Start by reading in information on the AO basis set as well as atomic positions
  OPEN(35,file='basis_set.info')

  READ(35,'(A)')comment
  WRITE(6,*)comment

  READ(35,*)

  READ(35,*)n_atom
  READ(35,*)s_dim

  WRITE(6,*)'n_atom',n_atom,'s_dim',s_dim

  READ(35,*)
  DO j=1,3
     READ(35,*)a(:,j)
  ENDDO
  !a = ABS(a)
  WRITE(6,*)'lattice vectors in bohr'
  WRITE(6,*)a


  ALLOCATE(atom_pos(n_atom,3),atomic_number(n_atom))

  READ(35,*)

  DO k=1,n_atom
     READ(35,*)atomic_number(k)
     READ(35,*)atom_pos(k,:)
  ENDDO

  READ(35,*)

  ALLOCATE(AO_basis(s_dim))


  DO nu=1,s_dim
     READ(35,*)AO_basis(nu)%num_gauss
     ALLOCATE(AO_basis(nu)%alpha(AO_basis(nu)%num_gauss),AO_basis(nu)%coeff(AO_basis(nu)%num_gauss),AO_basis(nu)%norm(AO_basis(nu)%num_gauss))
     READ(35,*)AO_basis(nu)%alpha
     READ(35,*)AO_basis(nu)%coeff
     READ(35,*)AO_basis(nu)%norm
     READ(35,*)AO_basis(nu)%pos
     READ(35,*)AO_basis(nu)%ncart
     ALLOCATE(AO_basis(nu)%cart_coeff(AO_basis(nu)%ncart),AO_basis(nu)%cart_mat(AO_basis(nu)%ncart,3))
     DO k=1,AO_basis(nu)%ncart
        READ(35,*)AO_basis(nu)%cart_coeff(k),AO_basis(nu)%cart_mat(k,:)
     ENDDO
     READ(35,*)
  ENDDO
  CLOSE(35)

  !DO nu=1,s_dim
  !   WRITE(6,*)'basis function',nu
  !   WRITE(6,*)AO_basis(nu)%num_gauss
  !   WRITE(6,*)AO_basis(nu)%alpha
  !   WRITE(6,*)AO_basis(nu)%coeff
  !   WRITE(6,*)AO_basis(nu)%norm
  !   WRITE(6,*)AO_basis(nu)%pos
  !   WRITE(6,*)AO_basis(nu)%ncart
  !   DO k=1,AO_basis(nu)%ncart
  !      WRITE(6,*)AO_basis(nu)%cart_coeff(k),AO_basis(nu)%cart_mat(k,:)
  !   ENDDO
  !   WRITE(6,*)
  !ENDDO


  !Now read in the coefficients of the NBOs in the AO basis described above
  !This file also contains information on the indices of the unit cells used in the calculation

!!!TRG Allowing different input files for different types of bonds; 
!!!TRG the filename to be specified as the first command line argument
 
!  OPEN(45,file='nbo_vis.out')
  call getarg (1, inputfile)
  open (45,file=inputfile) 
 
  READ(45,'(A)')comment
  WRITE(6,*)comment

  READ(45,*)

  READ(45,*)ao_nbasis
  WRITE(6,*)'ao_nbasis',ao_nbasis

  IF(ao_nbasis /= s_dim)THEN
     WRITE(6,*)'the number of basis functions from basis_set.info and nbo_vis.out dont match, are they for the same system?'
     STOP
  ENDIF

  READ(45,*)nnbo
  READ(45,*)num_l
  READ(45,*)nspins

  WRITE(6,*)'nnbo',nnbo,'num_l',num_l

  ALLOCATE(ao_coeff(s_dim,nnbo,num_l,nspins),index_l(num_l,3))

  READ(45,*)

  DO il=1,num_l
     READ(45,*)index_l(il,:)
  ENDDO

  READ(45,*)

  DO ispin=1,nspins
  DO il=1,num_l
     DO inbo=1,nnbo
        READ(45,*)ao_coeff(:,inbo,il,ispin)
     ENDDO
     READ(45,*)
  ENDDO
  ENDDO

  CLOSE(45)


  !The density is potentaily desired outside of the central unit cell, for instance a bind that spans unit cells
  !Thus, we will numerically integrate a grid over the volume 'box' which will contain the central unit cell and surrounding space.
  !The value of box is the absolute value of the lengths of the box
  !'origin' then account for the shift of the grid origin off center

  !This is the default box for surrounding space around the central unit cell
!  box = 3.d0*a
!  origin = -1.5d0*(a(:,1)+a(:,2)+a(:,3))
  !This is the default box for the central unit cell.  This can be used for supercells, where only the center of the box is of interest
!  box = a
!  origin = 0.d0

  !The mexh parameter determines how many grid points are used in each direction for the numerical integration
  !For non-cubic boxes, different components of this vector should be assigned so that resoultion is similar in all directions.
!  mesh = 20
 ! mesh (1) = 180
  !mesh(3) = 200



!!!TRG allowing the box type and the mesh be specified by the user through optional
!!!TRG command line arguments
if (iargc()>1) then
  call getarg(2, arg2)
  call getarg(3, arg3)
  call getarg(4, arg4)
  read (arg2,*) mesh(1)
  read (arg3,*) mesh(2)
  read (arg4,*) mesh(3)
  if (iargc()>4) then
    call getarg(5,arg5)
    read (arg5,*) ifperiodic
  else 
    ifperiodic=.true.
  end if
else 
  mesh=50
  ifperiodic=.true.
end if

if (ifperiodic) then
box=3.d0*a
origin=-1.5d0*(a(:,1)+a(:,2)+a(:,3))
else
box=a
origin=0.d0
end if



 

  WRITE(6,*)'box'
  WRITE(6,*)box
  WRITE(6,*)'origin',origin

  !ALLOCATE(density(nnbo),wave_funct(nnbo))
  ALLOCATE(norm_test(nnbo,nspins))
  norm_test = 0.d0

  !Prepare the file names of all the .cube files for each NBO of each type of spin
  !If there is only one spin type, it is unecessary to give different names based on spoin type

!!!TRG changing the output file naming:
  ALLOCATE(file_names(nnbo,nspins))
!TRG  IF( nspins .EQ. 1 )THEN
     DO inbo=1,nnbo
        IF( inbo < 10 )THEN
               write(file_names(inbo,1),'(A1, I1, A5)') '_',inbo,'.cube'
               file_names(inbo,1)=trim(inputfile(5:(len_trim(inputfile)-5)))//trim(file_names(inbo,1))
!TRG           WRITE(file_names(inbo,1),'(A5,I1,A5)')'nbo_',inbo,'.cube'
        ELSEIF( inbo < 100 )THEN
               write(file_names(inbo,1),'(A1, I2, A5)') '_',inbo,'.cube'
               file_names(inbo,1)=trim(inputfile(5:(len_trim(inputfile)-5)))//trim(file_names(inbo,1))
!
!TRG           WRITE(file_names(inbo,1),'(A5,I2,A5)')'nbo_',inbo,'.cube'
        ELSEIF( inbo < 1000 )THEN
               write(file_names(inbo,1),'(A1, I3, A5)') '_',inbo,'.cube'
               file_names(inbo,1)=trim(inputfile(5:(len_trim(inputfile)-5)))//trim(file_names(inbo,1))
!TRG           WRITE(file_names(inbo,1),'(A5,I3,A5)')'nbo_',inbo,'.cube'
        ELSE
           WRITE(6,*)'the code is not set up to write out the density for more than 999 bands'
           STOP
        ENDIF
     ENDDO
!TRG  ELSE  !For spin polarized calculations call the first spin type alpha and the second type beta NBO's

!TRG     DO inbo=1,nnbo
!TRG        IF( inbo < 10 )THEN
!TRG           WRITE(file_names(inbo,1),'(A11,I1,A5)')'alpha_nbo_',inbo,'.cube'
!TRG           WRITE(file_names(inbo,2),'(A10,I1,A5)')'beta_nbo_',inbo,'.cube'
!TRG        ELSEIF( inbo < 100 )THEN
!TRG           WRITE(file_names(inbo,1),'(A11,I2,A5)')'alpha_nbo_',inbo,'.cube'
!TRG           WRITE(file_names(inbo,2),'(A10,I2,A5)')'beta_nbo_',inbo,'.cube'
!TRG        ELSEIF( inbo < 1000 )THEN
!TRG           WRITE(file_names(inbo,1),'(A11,I3,A5)')'alpha_nbo_',inbo,'.cube'
!TRG           WRITE(file_names(inbo,2),'(A10,I3,A5)')'beta_nbo_',inbo,'.cube'
!TRG        ELSE
!TRG           WRITE(6,*)'the code is not set up to write out the density for more than 999 bands'
!TRG           STOP
!TRG        ENDIF
!TRG     ENDDO
!TRG  ENDIF

  !Calculates the grid volume for the specified mesh dimensions
  !This is used in numerically integrating the density
  CALL grid_volume(gridvol, box, mesh)
  WRITE(6,*)'integration grid volume',gridvol
  WRITE(6,*)

  !To aid in visualization we will want more than just the atoms in the unit cell.
  !To do this we will scan over all atom positions in all unit cells, and store those 'close' to the central unit cell
  !Start by establishing cutoffs for determining 'closeness'
  image_screen = 0.d0
  DO j=1,3
    image_screen(j,1) = MIN( origin(j), MINVAL(0.5d0*box(j,:)))
    image_screen(j,2) = MAXVAL(box(j,:)) + origin(j)
  ENDDO

  WRITE(6,*)'min image_screen',image_screen(:,1)
  WRITE(6,*)'max image_screen',image_screen(:,2)

  !Then every position is scanned over (each atom in all unit cells characterized by an l-vector)
  !Those within the cutoffs determined above are stored
  ALLOCATE(image_disp(n_atom*num_l))
  image_disp = .FALSE.
  image_count = 0
  image_tot = 0
  DO j=1,n_atom
     !WRITE(6,*)'atom',j,atom_pos(j,:)
     DO il=1,num_l
        image_tot = image_tot + 1
        image_pos = atom_pos(j,:)
        DO k=1,3
           image_pos = image_pos + index_l(il,k)*a(:,k)
        ENDDO
        IF( image_pos(1) > image_screen(1,1) .AND. image_pos(1) < image_screen(1,2))THEN
          IF( image_pos(2) > image_screen(2,1) .AND. image_pos(2) < image_screen(2,2))THEN
            IF( image_pos(3) > image_screen(3,1) .AND. image_pos(3) < image_screen(3,2))THEN
              WRITE(6,*)index_l(il,:)
              WRITE(6,*)image_pos
              image_disp(image_tot) = .TRUE.
              image_count = image_count + 1
            ENDIF
          ENDIF
        ENDIF
     ENDDO
     WRITE(6,*)
  ENDDO
  WRITE(6,*)'image count',image_count
  WRITE(6,*)

  !Then a file containing the dimensions of the central unit cell as well as the positions of all atoms in the central and surrounding unit cells
  OPEN(10, file='lattice_vec.cube')
  WRITE(10, *) 'Cube file generated by write_cube subroutine'
  WRITE(10, *) 'Density'
  WRITE(10, '(I5,3F12.6)' ) image_count, 0.d0, 0.d0, 0.d0
  WRITE(10, '(I5,3F12.6)' ) mesh(1), a(1,1)/mesh(1), a(2,1)/mesh(1), a(3,1)/mesh(1)
  WRITE(10, '(I5,3F12.6)' ) mesh(2), a(1,2)/mesh(2), a(2,2)/mesh(2), a(3,2)/mesh(2)
  WRITE(10, '(I5,3F12.6)' ) mesh(3), a(1,3)/mesh(3), a(2,3)/mesh(3), a(3,3)/mesh(3)
  image_tot = 0
  DO j=1,n_atom
     DO il=1,num_l
        image_tot = image_tot + 1
        IF( image_disp(image_tot) )THEN
           image_pos = atom_pos(j,:)
           DO k=1,3
              image_pos = image_pos + index_l(il,k)*a(:,k)
           ENDDO
           WRITE(10, '(I5,4F12.6)' ) atomic_number(j), 1.0, image_pos(1), image_pos(2), image_pos(3)
        ENDIF
     ENDDO 
  ENDDO

  !The remainder of the lattice_vec.cube file is filled with appropriate dimension and formatted zero-valued density
  !This prevents VMD from having trouble opening the file and displaying the appropriate unit cell vectors
  DO j=1,PRODUCT(DBLE(mesh))
     WRITE(10,"(E13.5)",ADVANCE='NO')0.d0
     IF( MODULO(MODULO(j,mesh(3)),6) == 0 )WRITE(10,*)
  ENDDO

  CLOSE(10)

  !Finally a .cube file for each nbo, for each spin is filled with appropriate density on the 'mesh' grid
  !Note: this grid can extend beyond the unit cell based on the box variable


  DO ispin=1,  nspins
!!!TRG     WRITE(6,*)'For spin type',ispin

     DO inbo=1,nnbo
        WRITE(6,*)'orbital number',inbo

        !Start by opening the appropriate file (names created above) and writing the cube file header 
        file_num = inbo+6
        OPEN(file_num, file=file_names(inbo,ispin))
        write(file_num, *) 'Cube file generated by write_cube subroutine'
        write(file_num, *) 'Density'
        write(file_num, '(I5,3F12.6)' ) n_atom, origin(1), origin(2), origin(3)
        write(file_num, '(I5,3F12.6)' ) mesh(1), box(1,1)/mesh(1), box(2,1)/mesh(1), box(3,1)/mesh(1)
        write(file_num, '(I5,3F12.6)' ) mesh(2), box(1,2)/mesh(2), box(2,2)/mesh(2), box(3,2)/mesh(2)
        write(file_num, '(I5,3F12.6)' ) mesh(3), box(1,3)/mesh(3), box(2,3)/mesh(3), box(3,3)/mesh(3)
        DO j=1,n_atom
           WRITE(file_num, '(I5,4F12.6)' ) atomic_number(j), 1.0, atom_pos(j,1), atom_pos(j,2), atom_pos(j,3)
        ENDDO

        !Then the grid is looped over for the total box and density tabulated
        DO ix=1,mesh(1)
           DO iy=1,mesh(2)
              DO iz=1,mesh(3)
                 density = 0.d0

                 r = ((ix-1)/DBLE(mesh(1)))*box(:,1) + ((iy-1)/DBLE(mesh(2)))*box(:,2) + ((iz-1)/DBLE(mesh(3)))*box(:,3)
                 r = r + origin

                 wave_funct = 0.d0
                 !The effect of each basis function in each unit cell must be calculated for each grid point
                 DO nu=1,s_dim
                    DO il=1,num_l

                       !In many unit cells, the coefficient is zero, so there is no need to calculate a function's value
                       IF( ao_coeff(nu,inbo,il,ispin) /= 0.d0 )THEN
                          r_pos = r - (AO_basis(nu)%pos + index_l(il,1)*a(:,1) + index_l(il,2)*a(:,2) + index_l(il,3)*a(:,3))
                          rsqr = DOT_PRODUCT(r_pos,r_pos)

                          !The guassian component of nu at r is now calculated. Screening is used
                          gaussian = 0.d0
                          DO j=AO_basis(nu)%num_gauss,1,-1  !In going backwards I assume alpha's are stored in descending order
                             screen = AO_basis(nu)%alpha(j) * rsqr
                             IF( screen .GT. 45.d0 )EXIT    !In the loop each gaussian is less diffuse than the last, so all following gaussians will also be screened
                             gaussian = gaussian + AO_basis(nu)%coeff(j)*AO_basis(nu)%norm(j)*EXP(-screen)
                          ENDDO

                          !The cartesian component and multiplication with the coefficient of the AO in the NBO is only done, it the gaussian component is non-zero
                          !Since many gaussians are screened out because of length, this screens alot of unnecessary calculations
                          IF( gaussian .NE. 0.d0 )THEN
                              cartesian = 0.d0
                              DO j=1,AO_basis(nu)%ncart
                                 cartesian = cartesian + AO_basis(nu)%cart_coeff(j)*PRODUCT(r_pos**AO_basis(nu)%cart_mat(j,:))
                              ENDDO
                              wave_funct = wave_funct + gaussian*cartesian*ao_coeff(nu,inbo,il,ispin)
                          ENDIF
                       ENDIF

                    ENDDO
                 ENDDO

                 !The wave function value is then squared to obtain a density
                 IF( wave_funct .NE. 0.d0 )THEN
                    density = wave_funct!**2
                    norm_test(inbo,ispin) = norm_test(inbo,ispin) + density*gridvol
                    IF( ABS(density) .LT. 1.d-30 )density = 0.d0  !Too low of values in the cube file mess up VMD.  Since they are small anyways I round down to zero
                 ENDIF

                 !The density is then written in the appropriate format of a .cube file
                 WRITE(file_num, "(E13.5)", advance="no") density
                 IF ( MOD(iz,6) == 0 )THEN
                    WRITE(file_num,*)
                 ENDIF

              ENDDO  !End of loop in mesh through fastest direction (z)

              WRITE(file_num,*)

           ENDDO
        ENDDO  !End of loop over mesh

        CLOSE(file_num)

        WRITE(6,*)'norm test for nbo',inbo
        WRITE(6,*)SNGL(norm_test(inbo,ispin))
        WRITE(6,*)

     ENDDO  !End of loop over band

  ENDDO


ENDPROGRAM visual


FUNCTION cross(a, b)
  REAL*8,INTENT (in) :: a(3), b(3)
  REAL*8 :: cross(3)

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross

SUBROUTINE grid_volume(gridvol, b, mesh)
  IMPLICIT NONE
  REAL*8, INTENT(IN)    ::   b(3,3)          ! Reciprocal lattice vectors
  INTEGER, INTENT(IN)   ::   mesh(3)         ! Number of grid points in each direction, inverse gives lenght of grid vector
  REAL*8, INTENT(OUT)   ::   gridvol         ! Volume of box to use for integration
  REAL*8                ::   grid_vec(3,3)    ! Vectors of gridbox in direction of recip lattice vectors
  INTEGER               ::   j

  INTERFACE
     FUNCTION cross(a,b)
       REAL*8,INTENT (in) :: a(3), b(3)
       REAL*8 :: cross(3)
     END FUNCTION cross
  END INTERFACE


  DO j=1,3
  grid_vec(:,j)=(1.d0/DBLE(mesh(j)))*b(:,j)
  ENDDO

  gridvol = DOT_PRODUCT(grid_vec(:,1), cross(grid_vec(:,2),grid_vec(:,3)))


END SUBROUTINE grid_volume




