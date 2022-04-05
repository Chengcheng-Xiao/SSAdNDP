
MODULE ssadndp_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: SSAdNDP

  CONTAINS
  
    SUBROUTINE SSAdNDP (inp, nspin, inputfile)
      USE nbo_shared
      USE matutil
      USE periodic_matutil
      
      IMPLICIT NONE
      type (nbo_input) :: inp
      integer :: nspin, ispin  
      character (15) :: inputfile
      
      ! write files for visualization?
      logical :: ifvis

      ! unit cell volume from CONTCAR
      real*8 :: vol
      ! lattice vectors and scaled
      real*8, dimension (3,3) :: lvecs, scaledlvecs
      ! atom coordinates from CONTCAR
      real*8, dimension (inp%natom,3) :: atcoords 
      ! determines the type of coordinates used in CONTCAR
      logical :: ifdirect
      ! number of atoms per unit cell (=inp%natom), number of atoms in the 'central'+ 26 (for 3D systems) neighbouring unit cells
      integer :: nat, nat27
      ! remapping indices from inp%indexg, preliminarily modified in projection? to include the 'central' and two neighbouring cells in every direction, back to the one-neighbour-system only.
      ! number of unit cells in the one-neighbour-system
      integer :: newng
      ! indices (in index%ng) of the neighbouring cells only
      integer, dimension (:), allocatable :: neighbournums
      ! counters, iterations
      integer :: i,j
      ! number of bonds for a given n
 

      ! general search:
        ! if read/write residual density (RD) file
        logical :: ifinputok !, ifreadRD, ifwriteRD
        ! number of different types of bonds in general search 
        !(different "n"s to try for nc-2e bonds) 
        integer ::  iset
        integer, dimension (nspin) :: nsets
        ! number of centers, ON thresholds, Distance thresholds
        integer, dimension (:, :), allocatable :: nctrs
        real*8, dimension (:, :), allocatable :: ONthrs, Dthrs
        !.true. to write all combns to combinations.out
        logical, parameter :: writecombns=.false. 
        !.true. to write all combinations of atomic centers (with at least one
        ! atom from the 'central' unit cell) and distance tests to disttest.out
        logical, parameter :: writedistres=.false. 
        !.true. to write to screen all comb-nations that passed testdist and testperiod
        ! and are submitted to blcomb and so on.
        logical, parameter :: writetestedcntrs=.false.
        integer :: cntbond
        ! maximum number of bonds possible for any n; presented as a parameter 
        ! since the bonds are found 'on-the fly' in gen_search; same for coeffs. 
        integer, parameter :: maxnbonds=1000
        integer, parameter :: maxcoeffs=3000
        ! Depletedensity uses actual sizes for nctrbonds and nctrONs;
        ! but dim 1 of nctrcoeffs will use the size of the largest subblock for an iset.
        integer:: maxsubblocksize
        ! the revealed bonds for a given n
        integer, dimension (:,:,:), allocatable :: nctrbonds
        ! their ONs:
        real*8, dimension (:), allocatable :: nctrONs 
        ! and coeffs:
        real*8, dimension (maxcoeffs, maxnbonds) :: nctrcoeffs       
        ! combination of nctrs(iset, nspin) numerals in range 1..nat27  
        integer, dimension (:), allocatable :: combn

       ! user-directed search:
        ! number of centers in UD
        integer, dimension (:), allocatable :: UDnctrs
        ! specification of centers for UD-search
        integer, dimension (:,:,:), allocatable :: UDctrs
        ! number of bonds
        integer :: UDnbonds      
        ! fragments, bonds, and centers iterators
        integer :: ifr, ibond, ictr
        ! density mtrx subblock
        real*8, dimension (:,:), allocatable :: UDsubblock
        integer :: UDsubblocksize
        ! occupations
        real*8, dimension (:), allocatable :: UDoccns
        ! coefficients
        real*8, dimension (:,:), allocatable :: UDcoeffs
        ! track first atom number
        integer :: firstat



      nat=inp%natom
 
!****** remapping indices from 'expanded' inp%indexg to neighbours only  

!!!check this with the dowloaded projection; whether indexg (inp%ng) has 25 unit cells (j's) -planar with two rows
!write (*, *) 'inp%indexg'
!write (*, *)
!write (*, *) 'inp%ng', inp%ng
!write (*,'(25I3)') (i, i=1,inp%ng)
!write (*, *) 
!write (*,'(25I3)') ((inp%indexg(i,j), j=1,inp%ng), i=1,3)

      j=0
      do i=1,3
          if (maxval (inp%indexg (i,:))>0) j=j+1
      end do
      newng=3**j
!write (*,*) 'newng', newng
      allocate (neighbournums(newng))
      j=0

      do i=1, inp%ng
        if ((abs(inp%indexg(1,i))<2).and.(abs(inp%indexg(2,i))<2).and.(abs(inp%indexg(3,i))<2))&
            then
            j=j+1
            neighbournums(j)=i
        end if
      end do

      nat27=nat*newng
!write (*,*) 'nat27', nat27

      call read_coords (vol, lvecs, atcoords, nat, ifdirect)
!write (*,*), vol, lvecs, atcoords, nat, ifdirect
      if (.not. ifdirect) then
        write (*,*) 'coordinates in CONTCAR are not Direct? Stopping'
        STOP
      end if
      scaledlvecs=lvecs*vol

      if (inputfile=='SSAdNDP.inp') then
        write (*,*) 'Starting general SSAdNDP search'
        write (*,*)
        call gen_search
      else if (inputfile=='SSAdNDP-UD.inp') then
        write (*,*) 'Starting user-directed SSAdNDP search'
        write (*,*)
        call ud_search
      end if



      CONTAINS

        SUBROUTINE gen_search

          IMPLICIT NONE
          integer, dimension (:), allocatable :: temp
          !visualization filename

          character (30) :: visfilename
          character (20) :: tempname

          open (54, file='SSAdNDP.out')
          write (54,*) 'General SSAdNDP Search'
          write (54,*)

!         call read_input (ifreadRD, ifwriteRD, nsets, nctrs, ONthrs, Dthrs, ifinputok)


!          if (ifreadRD) then

!            open (64, file='residmtrx.dat', form="unformatted")
!            read (64) inp%rho0
!            close (64)
!          end if
    

          do ispin=1, nspin
!           if (ifreadRD) then
!             if (nspin==2) then
!               write (54,'(A47, I1, A3, F15.10)')&
!                     'Initial density (from residmtrx.dat) for spin',&
!                     ispin,': ', resid()
!             else
!               write (54,*) 'Initial density (from residmtrx.dat):', resid()
!             end if
!           else
              if (nspin==2) then
                write (54,'(A26, I1, A2, F15.10)') ' Initial density for spin ',&
                       ispin,': ', resid()
              else
                write (54,*) 'Initial density:', resid()
              end if
!           end if
          end do

          call read_input (nsets, nctrs, ONthrs, Dthrs, ifinputok)
          if (.not. ifinputok) then
            write (*,*) 'ERROR: smth must be wrong with the input file'
            write (*,*) 'Check SSAdNDP.inp'
            write (*,*) 'Stopping'
            STOP
          end if

      
          if (writecombns) open (74, file='combinations.out')
          if (writedistres) open (75, file='disttest.out')
      

          do ispin=1, nspin
            write (54,*)
            write (54,*)

            if (nspin==2) write (54,'(A15, I1, A9)') ' ******** SPIN ', &
                                ispin, ' ********'

            write (54,*)  
            write (54,*) '----------------------------------------------------------'
            if (nspin==2) then
               write (54,'(A16, I1)') 'INPUT for spin ', ispin
            else 
               write (54,*) 'INPUT:'
            end if

            write (54,*)
            write (54,*) 'Number of different types of bonds in general search:'
            write (54, *)  nsets(ispin)  

            write (54,*) 'Number of centers, ON threshold, Distance threshold:'
            do iset=1, nsets(ispin)
              !write (*,*) nctrs(iset, ispin), nspin
              write (54,*) nctrs(iset, ispin), ONthrs(iset, ispin), Dthrs(iset,ispin)
            end do
            write (54,*) '----------------------------------------------------------'
            write (54,*)

            do iset=1, nsets(ispin) 

              write (*,*)          
              write (*,'(I3, A19)') nctrs(iset,ispin), '-center bond search'    
              write (*,*)       

              cntbond=0
              maxsubblocksize=0
              allocate (combn(nctrs(iset, ispin))) 
              if (writecombns) write (74, *) 'combinations of', nctrs(iset, ispin), &
                 'numerals out of ', nat27
              allocate (nctrbonds(maxnbonds, 2, nctrs(iset, ispin)))
              allocate (nctrONs(maxnbonds))

              call bond_search (1, nctrs(iset, ispin), Dthrs(iset, ispin), ONthrs(iset, ispin), &
                                cntbond, nctrbonds, nctrONs, nctrcoeffs, maxsubblocksize) 
       
              if (cntbond.ne.0) then
                write (54,*)
                write (54,'(I3, A14)') nctrs(iset,ispin), '-center bonds:' 
                call sortONs (nctrONs(1:cntbond), nctrbonds(1:cntbond,1:2,1:nctrs(iset, ispin)), &
                             nctrcoeffs(1:maxsubblocksize,1:cntbond)) 

                
                do ibond=1, cntbond
                  write (54,'(I6, A9, F20.15 )') ibond, ')     ON=',&
                         nctrONs(ibond)
                  write (54,*), '    atoms  (cells)'
                  write (54, 14) (nctrbonds(ibond,1,i), '      (', &
                        inp%indexg(:, nctrbonds(ibond,2,i)), ')' , i=1, nctrs(iset,ispin))
                  14 FORMAT (I6, A7, 3I3, A1)
                  write (54, *)
                end do
 

                call depletedensity (nctrbonds(1:cntbond,1:2,1:nctrs(iset, ispin)), cntbond, &
                                     nctrcoeffs(1:maxsubblocksize,1:cntbond), &
                                     nctrs(iset, ispin), nctrONs(1:cntbond))
                write (54, *) 'residual density:', resid()

              
                if (ifvis) then
                  allocate (temp(cntbond))
                  ! this is not really necessary here, but was introduced so that the same
                  ! visualization procedure can be used for both general and UD searches
                  temp=nctrs(iset, ispin)
                  ! name vis.out file
                  
                  write (tempname, '(I4,A11)') nctrs(iset, ispin), 'c-bonds.out'
                  tempname=adjustl(tempname)
                
                  if (nspin==2) then
                    if (ispin==1) visfilename= 'vis_gs_spin-1_'//tempname
                    if (ispin==2) visfilename= 'vis_gs_spin-2_'//tempname            
                  else
                    visfilename='vis_gs_' //tempname
                  end if
                  visfilename=adjustl(visfilename)
                  write (*,*) 'Preparing ', visfilename
                  write (*,*)
                  open (95, file=visfilename)


                  call visualization (nctrbonds(1:cntbond,1:2,1:nctrs(iset, ispin)), &
                                    cntbond, temp,nctrcoeffs(1:maxsubblocksize,1:cntbond))
                  deallocate (temp)
                  close (95)    
                end if

              else
                write (54,*)
                write (54,'(I3, A25)') nctrs(iset,ispin), '-center bonds: none found'
              end if

              deallocate (combn)
              deallocate (nctrbonds)
              deallocate (nctrONs)
            end do      
          end do

          if (writecombns) close (74)
          if (writedistres) close (84)
      
          close (54)
          

 !        if (ifwriteRD) then
 !    
 !          open (64, file='residmtrx.dat', form="unformatted", position="rewind")
 !          write (64) inp%rho0
 !          close (64)
   
 !        end if

          deallocate (nctrs)
          deallocate (ONthrs)
          deallocate (Dthrs)
    
        END SUBROUTINE gen_search

        SUBROUTINE sortONs (ons, bonds, cfs)

        IMPLICIT NONE
        real*8, dimension (:), intent(inout) :: ons
        real*8, dimension (size(ons)) :: tempons
        real*8, dimension (:,:), intent(inout) :: cfs
        real*8, dimension (size(cfs,1),size(cfs,2)):: tempcfs
        integer, dimension (:,:,:), intent(inout) :: bonds
        integer, dimension (size(bonds,1),size(bonds,2),size(bonds,3)) :: tempbonds
        integer, dimension (1) :: oldindex
        integer :: i
        real*8 :: largest
 
        tempons=ons
        tempcfs=cfs
        tempbonds=bonds
        cfs=0
        bonds=0
        ons=0

        
        do i=1, size (ons) 

        oldindex=MAXLOC(tempons)
        ons(i)=tempons(oldindex(1))
        tempons(oldindex(1))=0.0
        bonds(i,:,:)=tempbonds(oldindex(1),:,:)
        cfs(:,i)=tempcfs(:,oldindex(1))
!       write (54,*) oldindex(1)

        end do
        write (54,*) 
        
        
        END SUBROUTINE sortONs

        SUBROUTINE ud_search

        IMPLICIT NONE
        real*8, dimension (:,:), allocatable :: UDallcoeffs
        integer :: maxUDsubblocksize, locn, f, istat
        character (30) :: visfilename
        !counter of bonds 
        integer :: iUDbond

        maxUDsubblocksize=0

        open (17, file='SSAdNDP-UD.inp', status='old', iostat=istat)

        if (istat.ne.0) then
          write (*,*) 'could not open input file, SSAdNDP-UD.inp'
        end if
        read (17,*)

        read (17,*)
        read (17,*, iostat=istat) ifvis
        if (istat.ne.0) then
          write (*,*) 'ERROR: smth must be wrong with the input file'
          write (*,*) 'Check SSAdNDP-UD.inp'
          write (*,*) 'Stopping'
          STOP
        end if

        open (55, file='SSAdNDP-UD.out')
        write (55,*) 'User-Directed Search:'
        write (55,*) 


        do ispin=1, nspin
          if (nspin==1) then
              write (55, *) 'Initial residual density:', resid()
          else 
              write (55, '(A25, I1, A1, F15.10)') 'Initial density for spin ', ispin, ':', resid()
          end if
        end do

        do ispin=1, nspin
          write (55,*)


          call read_input_ud(UDnctrs, UDctrs, UDnbonds, ifinputok)

          if (.not. ifinputok) then
            write (*,*) 'ERROR: smth must be wrong with the input file'
            write (*,*) 'Check SSAdNDP-UD.inp'
            write (*,*) 'Stopping'
            STOP
          end if

          allocate (UDallcoeffs(maxcoeffs, UDnbonds))


          !do ifr=1, size (UDctrs, 1)
            !do ibond=1, UDnbonds(ifr)
         
          do iUDbond=1, UDnbonds

              write (55,*) ''
              write (55, 21) ' Fragment with ', UDnctrs(iUDbond), ' atomic centers'
              21 FORMAT (A15, I3, A15)
              write (55,*) ' atoms (cells):'
              write (55, 16) (UDctrs(iUDbond, 1, ictr), '  (', &
                             inp%indexg(:,UDctrs(iUDbond, 2, ictr)), ')', ictr=1,  UDnctrs(iUDbond))  
              16 FORMAT (I6, A3, 3I3, A1)
             
              call blcomb (UDnctrs(iUDbond), UDctrs(iUDbond, :, 1:UDnctrs(iUDbond)), &
                   UDsubblock, UDsubblocksize)  

              if (maxUDsubblocksize<UDsubblocksize) maxUDsubblocksize=UDsubblocksize
              allocate (UDoccns(UDsubblocksize))
              allocate (UDcoeffs(UDsubblocksize, UDsubblocksize))
          
              call matdiag(UDsubblock, UDoccns, UDcoeffs)

          
              write (55,*), 'ON=', maxval(UDoccns)
  
              do f=1, size (UDoccns,1)
                if (UDoccns(f)==maxval(UDoccns)) locn=f
              end do

              UDallcoeffs(1:UDsubblocksize,iUDbond)=UDcoeffs(1:UDsubblocksize, locn)
            
              call depletedensity (UDctrs(iUDbond, :, 1:UDnctrs(iUDbond)), 1, &
                                  UDcoeffs(:, maxloc(UDoccns)), &
                                  UDnctrs (iUDbond), UDoccns(maxloc(UDoccns)))
              write (55, *) 'residual density:', resid()

       
              deallocate (UDoccns)
              deallocate (UDcoeffs)          

!                write (tempname, '(I4,A11)') nctrs(iset, ispin), 'c-bonds.out'
!                tempname=adjustl(tempname)

          end do

 
          if (ifvis) then
     
            if (nspin==2) then
              if (ispin==1) visfilename= 'vis_ud_spin-1_bonds.out'
              if (ispin==2) visfilename= 'vis_ud_spin-2_bonds.out'          
            else
              visfilename='vis_ud_bonds.out'
            end if
            visfilename=adjustl(visfilename)

            open (95, file=visfilename)
            write (*,*) 'Preparing ', visfilename
            write (*,*) 


            call visualization (UDctrs(1:UDnbonds,1:2, 1:maxval(UDnctrs(1:UDnbonds))), UDnbonds, &
                             UDnctrs (1:UDnbonds), UDallcoeffs(1:maxUDsubblocksize, 1:UDnbonds))

            close (95)
          end if
  
          deallocate (UDnctrs)
          deallocate (UDctrs)
          deallocate (UDallcoeffs)
        end do

        close (17)


        close (55)

!       if (ifwriteRD) then
!         open (64, file='residmtrx.dat', form="unformatted", position="rewind")
!         write (64) inp%rho0
!         close (64)
!       end if

        END SUBROUTINE ud_search

        SUBROUTINE read_input_ud(UDnctrs, UDctrs, UDnbonds, ifinputok)


          IMPLICIT NONE

          logical, intent(out) :: ifinputok
          integer, dimension (:), allocatable, intent (out) :: UDnctrs
          integer, dimension (:,:,:), allocatable, intent (out) :: UDctrs
          integer, intent (out) :: UDnbonds
          integer :: x, y, z
          integer :: iUDctr, nfrs, ifr, ctr
          integer :: istat
          ! number of bonds on current fragment
          integer :: frnbonds, ifrbond
   
          UDnbonds=0

          ifinputok=.true.
          read (17,*)
          write (55, *) 
          if (nspin==2) write (55,'(A16, I1, A9)') ' ******** SPIN ', &
                                ispin, ' ********'
 
          write (55,*) 
          write (55,*) '----------------------------------------------------------'
          if (nspin==2) then
            write (55,'(A16, I1)') 'INPUT for spin ', ispin
          else
            write (55,*) 'INPUT:'
          end if
          write (55,*)
 
          read (17,*, iostat=istat) nfrs        
          write (55,'(A43,I3)') 'Number of different fragments to search on (spin 1):', nfrs
          if (istat.ne.0) then
            ifinputok=.false.
            goto 3344
          end if
          
          allocate (UDnctrs (maxnbonds))
          allocate (UDctrs(maxnbonds, 2, nat27))

          do ifr=1, nfrs

            UDnbonds=UDnbonds+1
            write (55,*) 
            read (17,*)
            read (17, *, iostat=istat) UDnctrs(UDnbonds)
            write (55, '(A30,I2,A2,I2)') ' Number of centers on fragment ', ifr, ': ',&
                                           UDnctrs(UDnbonds)
            if (istat.ne.0) then
              ifinputok=.false.
              goto 3344
            end if
            read (17,*)
            write (55, *) 'Atom numbers (cell a,b,c)'
            do iUDctr=1, UDnctrs(UDnbonds)
              read (17,*, iostat=istat) ctr, x, y, z
              write (55,'(I3, A3, 3I3, A1)') ctr, ' (', x, y, z, ')'
              if (istat.ne.0) then
                ifinputok=.false.
                goto 3344
              end if
              if ((ctr<1).or.(ctr>nat)) then
                write (*,*) 'smth wrong with Atom numbers in the input for user-directed search'
                ifinputok=.false.
                goto 3344
              end if
              if ((abs(x)>1).or.(abs(y)>1).or.(abs(z)>1)) then
                write (*,*) 'cell directions have to be within [-1,1]'
                ifinputok=.false.
                goto 3344
              end if
              UDctrs (UDnbonds, 1,iUDctr)=ctr
              UDctrs (UDnbonds, 2,iUDctr)=find_g(x,y,z)
            end do
            read (17,*)
            read (17,*, iostat=istat) frnbonds
            write (55, '(A14, I3, A20, I3)') 'Searching for', frnbonds, ' bond(s) on fragment', ifr
            if (istat.ne.0) then
              ifinputok=.false.   
              goto 3344
            end if
            if (frnbonds>1) then
              do ifrbond=1, frnbonds-1
                UDnbonds=UDnbonds+1
                UDnctrs(UDnbonds)=UDnctrs(UDnbonds-1)
                UDctrs (UDnbonds, 1:2, 1:UDnctrs(UDnbonds))=&
                         UDctrs (UDnbonds-1,1:2, 1:UDnctrs(UDnbonds-1))

              end do
            end if

          end do


          write (55,*) '----------------------------------------------------------'
          write (55,*)

3344    END SUBROUTINE read_input_ud

        SUBROUTINE visualization (nctrbonds, cntbond, cnctrs, nctrcoeffs)
          

          IMPLICIT NONE

          integer :: ibond, icell, ictr
          integer, intent (in) :: cntbond
          ! cnctrs is an array here, because when visualization is called, the
          ! number of centers is different in the UD search 
          integer, dimension (:), intent (in) :: cnctrs
          ! the revealed bonds for a given n
          integer, dimension (cntbond, 2, maxval(cnctrs)), intent(in) :: nctrbonds
          real*8, dimension (cntbond, (maxval(inp%ibasismap)-1 ),inp%ng) :: viscoeffs 
          integer :: kfirst, klast, ksize, kkfirst, kklast
          real*8, dimension (:, :), intent (in) :: nctrcoeffs       

  
          do ibond=1, cntbond
            viscoeffs(ibond, :,:)=0.d0
            kkfirst=1 

                
            do ictr=1, cnctrs (ibond)
              kfirst=inp%ibasismap(nctrbonds(ibond, 1, ictr))
              klast=inp%ibasismap(nctrbonds(ibond, 1, ictr)+1)-1
              ksize=klast-kfirst+1
              kklast=kkfirst+ksize-1

              do icell=1, inp%ng

                if (icell==nctrbonds(ibond, 2, ictr)) then
                  viscoeffs(ibond, kfirst:klast, icell)=nctrcoeffs(kkfirst:kklast, ibond)
                end if

              end do    
                 
              kkfirst=kklast+1  
                    
            end do
            
          end do
                 
          do ibond=1, cntbond
            viscoeffs(ibond,:,:) = periodic_matvecmul(inp%transform,viscoeffs(ibond,:,:))
         
          end do 

          !!! adapted from BDD's nbo_main.f90:

          write(95,*)"Output of lattice vector and ao coeffs from JRS periodic NBO code"
          write(95,*)
          write(95,*)inp%nbasis, '! number of basis functions per unit cell'  
          write(95,*)cntbond, '! number of possible lonepairs and NBOs, set of coefficients'
          write(95,*)inp%ng, '! number of l_vectors, unit cell pairs'
!         write(95,*)nspin,'! number of unique spins for which NBOs have been obtained'
          write(95,*)'1 ! number of unique spins is one, since spins 1 and 2 are in dif files'
          write(95,*)
   
          do icell=1,inp%ng
            write(95,*)inp%indexg(:,icell)
          end do

          write(95,*)

!          do ispin=1,nspin    
            do icell=1,inp%ng
              do ibond=1,cntbond
                write(95,*)viscoeffs(ibond, :, icell)
              end do
              write (95,*)
            end do
!          end do
           
     
        END SUBROUTINE visualization
    
    
        FUNCTION resid ()
       
          IMPLICIT NONE
          integer :: f
          real*8 :: resid
          resid=0.
          !write (*,*) nspin
          do f=1, size(inp%rho0, 1)
            resid=resid+inp%rho0(f,f,1,ispin)
          end do
        
        END FUNCTION resid


        SUBROUTINE depletedensity (nctrbonds, cntbond, nctrcoeffs, cnctrs, nctrONs)
       
!!!       USE matutil
!!!       USE periodic_matutil
          IMPLICIT NONE
          ! number of centers in current bond
          integer, intent (in) :: cnctrs
          ! number of bonds for a given n; 
          integer, intent (in) :: cntbond
          ! the revealed bonds for a given n
          integer, dimension (cntbond, 2, cnctrs), intent(in) :: nctrbonds
          ! their ONs:
          real*8, dimension (cntbond), intent(in) :: nctrONs  
          ! and coeffs:
          real*8, dimension (:, :), intent (in) :: nctrcoeffs       

          integer :: kfirst, klast, lfirst, llast, ksize, lsize, kkfirst, kklast, &
                     llfirst, lllast, kg, ibond,s, k ,l,t,x , sfirst, slast


          do ibond=1, cntbond

            kkfirst=1 

            do k=1, cnctrs

              kfirst=inp%ibasismap(nctrbonds(ibond, 1, k))
              klast=inp%ibasismap(nctrbonds(ibond, 1, k)+1)-1
              ksize=klast-kfirst+1
              kklast=kkfirst+ksize-1
              llfirst=1

              do l=1, cnctrs

                lfirst=inp%ibasismap(nctrbonds(ibond, 1, l))
                llast=inp%ibasismap(nctrbonds(ibond, 1, l)+1)-1
                lsize=llast-lfirst+1
                lllast=llfirst+lsize-1

                if (l==k) then

                  inp%rho0(kfirst:klast, lfirst:llast, 1, ispin)=&
                  inp%rho0(kfirst:klast, lfirst:llast, 1, ispin)-nctrONs(ibond)*&
                  outer_product(nctrcoeffs(llfirst:lllast,ibond),nctrcoeffs(llfirst:lllast,ibond))
                end if
            
                if (l>k) then

                  kg=find_g((inp%indexg(1,nctrbonds(ibond, 2, l))-&
                  inp%indexg(1,nctrbonds(ibond, 2, k))),(inp%indexg(2,nctrbonds(ibond, 2, l))-&
                  inp%indexg(2,nctrbonds(ibond, 2, k))),(inp%indexg(3,nctrbonds(ibond, 2, l))-&
                  inp%indexg(3,nctrbonds(ibond, 2, k))))
                    
            
                  if (kg==-1) then
                    write (*,*)'Error with indices in depletedensity'
                    STOP
                  end if 

                  inp%rho0(kfirst:klast, lfirst:llast, kg, ispin)=&
                  inp%rho0(kfirst:klast, lfirst:llast, kg, ispin)-nctrONs(ibond)*&
                  outer_product(nctrcoeffs(kkfirst:kklast,ibond),nctrcoeffs(llfirst:lllast,ibond))
                end if
  
                if (l<k) then
          
                  kg=find_g((inp%indexg(1,nctrbonds(ibond, 2, l))-&
                  inp%indexg(1,nctrbonds(ibond, 2, k))),(inp%indexg(2,nctrbonds(ibond, 2, l))-&
                  inp%indexg(2,nctrbonds(ibond, 2, k))),(inp%indexg(3,nctrbonds(ibond, 2, l))-&
                  inp%indexg(3,nctrbonds(ibond, 2, k))))
                   
                  if (kg==-1) then
                    write(*,*)'Error with indices in depletedensity'
                    STOP
                  end if

                  inp%rho0(kfirst:klast, lfirst:llast, kg, ispin)=&
                  inp%rho0(kfirst:klast, lfirst:llast, kg, ispin)-nctrONs(ibond)*&
                  outer_product(nctrcoeffs(kkfirst:kklast,ibond),nctrcoeffs(llfirst:lllast,ibond))
                end if
                  llfirst=lllast+1
              end do
                kkfirst=kklast+1  
            end do

!            write (54,'(A32, I3, A11, F15.10)')  '     Residual density after bond', ibond,&
!                  ') depleted:', resid ()

          end do
        
        END SUBROUTINE depletedensity  

        RECURSIVE SUBROUTINE bond_search (r, cnctrs, cDThr, cONThr, &
                             cntbond, nctrbonds, nctrONs, nctrcoeffs, maxsubblocksize)
  !       USE matutil
  !       USE nbo_shared
          IMPLICIT NONE
        
          ! for current set:
          ! number of atomic centers
          integer, intent (in) :: cnctrs 
          ! current thresholds:
          real*8, intent (in) :: cONThr, cDThr
          ! recursion
          integer :: r
          ! other counters
          integer :: inum, ictr, k
          ! number of bonds for a given n; 
          integer, intent (out) :: cntbond
          integer, dimension (2,cnctrs) :: atoms
          ! density mtrx subblock
          real*8, dimension (:,:), allocatable :: subblock
          integer :: subblocksize
          ! occupations
          real*8, dimension (:), allocatable :: occns
          ! coefficients
          real*8, dimension (:,:), allocatable :: coeffs
          ! the revealed bonds for a given n
          integer, dimension (maxnbonds, 2, cnctrs), intent(out) :: nctrbonds
          ! their ONs:
          real*8, dimension (maxnbonds), intent(out) :: nctrONs  
          ! and coeffs:
          real*8, dimension (maxcoeffs, maxnbonds), intent (out) :: nctrcoeffs
          integer, intent (out) :: maxsubblocksize       

        
          if (r > cnctrs) then    
          !! OUTER IF
            if (writecombns) write (74, *) combn

            !to include only comb-s with at least one atom from the central unit cell:
            if (combn(1)<=nat) then     
              if (combn(1).ne.firstat) then
                if (cnctrs>1) then
                  firstat=combn(1)
                  write (*,'(A33, I3)') 'testing combinations around atom ' , firstat
                  write (*,*)
                end if
              end if
            !! OUTER -1 IF
!             countaccept=countaccept+1
              atoms=assignatom(combn, cnctrs)
!             write(24,*), ''
!             write(24,*), 'atom   ',(wAts(1,k), k=1,NCent)
!             write(24,*), 'cell   ', (wAts(2,k), k=1,NCent)

              if ((testdist(atoms, cnctrs, cDThr)==.true.).and.(testperiod(combn,cnctrs)==.true.)) &
                then  
              !OUTER -2 IF
!!!             write (*,*)
                if (writetestedcntrs) then
                  write (*,*)  'testing centers:'
                  write (*,*), ' atoms    (cell locns)    [cell indx]:'
                  write (*, 114) (atoms(1,ictr), '     (', inp%indexg(:, atoms(2,ictr)), &
                                  ' )    [', atoms(2,ictr),']' , ictr=1, cnctrs)
                  114 FORMAT (I6, A6, 3I3, A7, I2, A1)
                  write (*,*)
                end if         

                
!!!             blcombmtrx=0.
!!!             blcombmtrxsize=0.
                call blcomb (cnctrs, atoms, subblock, subblocksize)

                allocate (occns(subblocksize))
                allocate (coeffs(subblocksize, subblocksize))
        
!!!             occ=0.
!!!             coeff=0.
                call matdiag (subblock, occns, coeffs)
                deallocate (subblock)
!!!             countaccept=countaccept+1

                do k=1, subblocksize
                  if (occns(k)>=cONThr) then
                  !OUTER -3 IF

                    write (*,*) 'accepting a bond:'
                    cntbond=cntbond+1


                    write (*, 13) cntbond,')', cnctrs, '-center bond. ON=', (occns(k))
                    13 FORMAT (I3, A1, I5, A17, F20.14)
                    write (*,*), ' atoms (cells):'
                    write (*, 14) (atoms(1,ictr), '  (', &
                    inp%indexg(:, atoms(2,ictr)), ')' , ictr=1, cnctrs)
                    14 FORMAT (I6, A3, 3I3, A1)
                    write (*, *)

                    if (cntbond>maxnbonds) then
                      write (*,*)
                      write (*,*) 'ERROR: Number of accepted bonds for the current nc-search'
                      write (*,'(1X, A32, I8)') 'is larger than the parametrized:', maxnbonds
                      write (*,*) 'If even a larger number is still reasonable,'
                      write (*,*) "change the 'maxnbonds' parameter"
                      write (*,*) 'in ssadndp_module.f90 and recompile.'
                      write (*,*) 'Otherwise, adjust the thresholds'
                      write (*,*) 'Now stopping'
                      stop
                    end if

                    nctrbonds(cntbond,1:2,1:cnctrs)=atoms(1:2,1:cnctrs)
                    nctrONs(cntbond)=occns(k)

                    nctrcoeffs(1:subblocksize,cntbond)=coeffs(1:subblocksize,k)
                    if (maxsubblocksize<subblocksize) maxsubblocksize=subblocksize
!                   write (*,*) 'maxsubblocksize', maxsubblocksize

                  end if
                  !!OUTER -3 ENDIF
                end do
                deallocate(occns)
                deallocate(coeffs)
              end if  
              !! OUTER -2 ENDIF
!           else    
            !! OUTER -1 ELSE
!             write(24,*), 'skipping this:', comb, '   combination'
            end if   
            !! OUTER -1 ENDIF
          else         
          ! OUTER ELSE 
            do inum = 1, nat27
              if ((r == 1) .or. (inum > combn (r - 1))) then
              !! OUTER -1A IF
                combn (r) = inum
!               print *, comb(m), '!comb(m)', 'where m=',m
                call bond_search (r+1, cnctrs, cDThr, cONThr, cntbond, nctrbonds, nctrONs, nctrcoeffs,maxsubblocksize)
              end if
              !! OUTER -1A ENDIF
            end do
          end if     
          ! OUTER END IF

        END SUBROUTINE bond_search


        SUBROUTINE blcomb (cnctrs, atoms, subblock, subblocksize)
!         USE matutil
!         USE periodic_matutil
 
          IMPLICIT NONE
          integer, intent (in) :: cnctrs 
          integer, dimension (2,cnctrs), intent (in) :: atoms
          integer, intent (out) :: subblocksize
          real*8, intent (out), dimension (:,:), allocatable :: subblock
          integer :: kfirst, klast, lfirst, llast, ksize,  kkfirst, kklast, llfirst, lllast, k, l, lsize, kg, x,y
          
!!!       kfirst=0
!!!       klast=0
!!!       lfirst=0
!!!       llast=0
!!!       kkfirst=0
!!!       kklast=0
!!!       llfirst=0
!!!       lllast=0

          subblocksize=0

          do k=1, cnctrs
            kfirst=inp%ibasismap(atoms(1,k))
            klast=inp%ibasismap(atoms(1,k)+1)-1
            ksize=klast-kfirst+1
            subblocksize=subblocksize+ksize
          end do
  
          allocate (subblock (subblocksize, subblocksize))
!!!       blcombmat=0.
    
          kkfirst=1
      
          do k=1, cnctrs
            kfirst=inp%ibasismap(atoms(1,k))
            klast=inp%ibasismap(atoms(1,k)+1)-1
            ksize=klast-kfirst+1
            kklast=kkfirst+ksize-1
            llfirst=1
    
            do l=1, cnctrs
              lfirst=inp%ibasismap(atoms(1,l))
              llast=inp%ibasismap(atoms(1,l)+1)-1
              lsize=llast-lfirst+1
              lllast=llfirst+lsize-1
            
              if (l==k) then

                subblock(kkfirst:kklast, llfirst:lllast)=inp%rho0(kfirst:klast, lfirst:llast, 1, ispin)
              end if
            
              if (l>k)  then

                kg=find_g((inp%indexg(1,atoms(2, l))-inp%indexg(1,atoms(2, k))),(inp%indexg(2,atoms(2, l))-&
                inp%indexg(2,atoms(2, k))),(inp%indexg(3,atoms(2, l))-inp%indexg(3,atoms(2, k))))
 
                if (kg==-1 ) then
                  write (*,*)'Error in looking for kg index in filling subblock. Stopping'
                  STOP
                end if  

                subblock (kkfirst:kklast, llfirst:lllast)=inp%rho0(kfirst:klast, lfirst:llast, kg, ispin)
              end if
  
              if (l<k)  then

                subblock(kkfirst:kklast,llfirst:lllast)=transpose(subblock (llfirst:lllast, kkfirst:kklast))
              end if

              llfirst=lllast+1
            end do
            kkfirst=kklast+1
          end do
     
        END SUBROUTINE blcomb


        FUNCTION distance (p1,p2)
     
          IMPLICIT NONE
          real*8, dimension (3), intent (in) :: p1, p2
          real*8, dimension (3) :: diff
          real*8 :: distance
      
          diff=p2-p1
          distance=sqrt(diff(1)**2+diff(2)**2+diff(3)**2)
          !write (*,*) 'distance', distance

        END FUNCTION distance


        FUNCTION cartesian (atomnum, cellnum)
      
          IMPLICIT NONE
          real*8, dimension (3) :: cartesian
          real*8, dimension (3,nat) ::  atcoords_tr
          integer, intent (in) :: atomnum, cellnum

          atcoords_tr=transpose(atcoords)
          cartesian=matmul(transpose(scaledlvecs), (atcoords_tr(:,atomnum)+inp%indexg(:,cellnum)))
          ! write (*,*) 'scaledlvecs', scaledlvecs

        END FUNCTION cartesian


        FUNCTION testperiod(atomnumslist, length)

          IMPLICIT NONE
          integer, intent (in) :: length
          integer :: x,y
          integer, dimension (length), intent (in) :: atomnumslist
          logical :: testperiod

          testperiod=.true.
          do x=1, length
            if (atomnumslist(x)> nat) then
              if (mod (atomnumslist(x), nat).ne.0) then
                if (mod(atomnumslist(x), nat)<atomnumslist(1)) then
                  testperiod=.false.
                  !write (*,*) 'testperiod false'
                end if
              end if
            end if
          end do

        END FUNCTION testperiod


!       FUNCTION go_positive(cells, length)

!!! commented out, since does not include all cells. substituted with testperiod.

!         IMPLICIT NONE
!         integer :: cells(:), w,z, length
!         logical :: go_positive

!         go_positive=.true.
!         do w=1,length
!           do z=1,3
!             if (inp%indexg(z,cells(w))<0) then
!               go_positive=.false.
!               exit
!             end if
!           end do
!         end do
 
!       END FUNCTION go_positive


        FUNCTION testdist (atomlist, length, thr)
 
          IMPLICIT NONE
          integer, intent (in) :: length
          integer, intent (in), dimension (2, length) :: atomlist
          integer :: s,r, nums(length)
          logical :: testdist
          real*8, dimension (3, length) :: crtsnpts 
          real*8, intent (in) :: thr
          real*8 ::  dist 
       
          testdist=.true.
          if (thr.ne.0) then  !!! IF-0
            do s=1, length
              crtsnpts(:,s)=cartesian(atomlist(1,s),atomlist(2,s))
            end do

            if (writedistres) then  !* IF-sep
              write(75,*) 'atom (cell)     cartesian coordinates:'
              write(75, 12) (atomlist (1,j),' (', inp%indexg(:,atomlist (2,j)),')', &
                            (crtsnpts(i,j), i=1,3), j=1,length)
              12 FORMAT (I3, A2, 3I2, A1, F14.2, 2F8.2)
            end if                  !* If-setp

            if (length > 1) then   !!! IF-1
              do s=1, length
                do r=s+1, length
                  dist=distance(crtsnpts(:,s),crtsnpts(:,r))

                  if (writedistres) then !* if-sep
                    write (75,*) 'distance b/w atom(cell) and atom(cell):'
                    write (75,11) atomlist (1,s), ' (', inp%indexg(:,atomlist (2,s)), ') and ', &
                                  atomlist (1,r), ' (', inp%indexg(:,atomlist (2,r)), '):', dist
                    11 FORMAT (I3, A2, 3I2, A6, I3, A2, 3I2, A2, F10.2)
                  end if        !* if-sep

                  if (dist>thr) then !!! IF-2
                    testdist=.false.
                    if (writedistres) then !*if-sep
                      write (75,*) 'combination not accepted' 
                      write (75,*) 
                    end if                 !*if-sep
                    goto 1234
                  end if             !!! IF-2
                end do
              end do 

            else                 !!! else-1
              if (writedistres) write (75,*) 'no distance check for 1 center'
            end if      !!! IF-1

           
            if (testdist==.true.) then
              if (writedistres) then
                write (75,*) 'combination accepted'
                write (75,*)
              end if
            end if

          else    !!! else-0 
            if (writedistres) then
              write (75,*) 'no distance check requested (input thr=0) for'
              write (75,'(I2, A13)') length, '-center bonds'
            end if
          end if   !!! if-0


1234    END FUNCTION testdist


        FUNCTION assignatom (numsarray, length)
       
          IMPLICIT NONE
          integer, intent (in) :: length
          integer, intent (in), dimension (length) ::numsarray
          integer :: t
          integer, dimension (2, length) :: assignatom     
     
          do t=1, length 
            assignatom(1,t)=mod(numsarray(t),nat)
            if (assignatom(1,t)==0) assignatom(1,t)=nat
          end do
          do t=1, length
            if (assignatom(1,t)==nat) then
              assignatom(2,t)=neighbournums(numsarray(t)/nat)
            else
              assignatom(2,t)=neighbournums(numsarray(t)/nat+1)
            end if
          end do

        END FUNCTION assignatom        



!       SUBROUTINE create_distmat(NAt, numcell, cellloc, distmat)

!!!     commented out since switched to 'on-the-fly' searching

!         IMPLICIT NONE
!         integer:: NAt, numcell, cellloc(3,numcell), m,n,o,p
!         real :: distmat(NAt,numcell,NAt,numcell), p1(3), p2(3)
!         p1=cartesian(1,1)
!         p2=cartesian(2,1)
!         print *, p1, 'p1'
!         print *, p2, 'p2'
!         print *, distance(p1,p2), 'distance p1-p2'
!         print *, distance(cartesian(1,1), cartesian(1,27)), '1,1 - 1,27'
!         print *, distance(cartesian(2,1), cartesian(2,27)), '2,1 - 2,27'
!         print *, distance(cartesian(2,1), cartesian(1,26)), '2,1 - 1,26'

!         do m=1, numcell
!           do n=1, NAt
!             do o=1, numcell
!               do p=1, NAt
!                 p1=cartesian(m,n)
!                 p2=cartesian(o,p)
!                 print *, p1!, p2
!               end do
!             end do
!           end do
!         end do

!       END SUBROUTINE create_distmat       


        SUBROUTINE read_input(nsets, nctrs, ONthrs, Dthrs, ifinputok)


          IMPLICIT NONE

          logical, intent (out) :: ifinputok !ifreadRD, ifwriteRD
          integer, dimension(nspin), intent (out) :: nsets
          integer, dimension (:, :), allocatable, intent (out) :: nctrs
          real*8, dimension (:, :), allocatable, intent(out) :: ONthrs, Dthrs
          integer :: iset, istat
          character*128 :: dataline

          write (*,*) 'Reading input file SSAdNDP.inp'
          write (*,*)          
          open (14, file='SSAdNDP.inp', status='old', iostat=istat)

          if (istat.ne.0) then
            write (*,*) 'could not open input file'
            ifinputok=.false.
            goto 2233
          end if

          read (14,*)
!         read (14,*)
          
!         read (14,*) ifreadRD

!         write (54,*) 
!         write (54,*) '----------------------------------------------------------'
!         write (54,*) 'INPUT:'
!         write (54,*)
!         write (54,*) 'Read from residual density file:', ifreadRD
!         read (14,*)
!         read (14,*) ifwriteRD
!         write (54,*) 'Update residual density file:', ifwriteRD
!         write (54,*)

          !using maxnbonds below, since nsets(ispin=1) and nsets(ispin=2) can be different
          allocate (nctrs(maxnbonds, nspin))
          allocate (ONthrs(maxnbonds, nspin))
          allocate (Dthrs(maxnbonds, nspin))
 
          ifinputok=.true.
       
          read (14,*)
          read (14,*, iostat=istat) ifvis
          if (istat.ne.0) then
            ifinputok=.false.
            goto 2233
          end if

          do ispin=1, nspin
            read (14,*) 
!           if (nspin==2) write (54,'(A6,I1)') ' Spin ', ispin
!           write (54,*) 'Number of different types of bonds in general search:'
            read (14, *, iostat=istat) nsets(ispin)
            if (istat.ne.0) then
              ifinputok=.false.
              goto 2233
            end if
!           write (54, *)  nsets(ispin)  

            read (14,*) 
!           write (54,*) 'Number of centers, ON threshold, Distance threshold:'
            do iset=1, nsets(ispin)
              read(14, *, iostat=istat) nctrs(iset, ispin), ONthrs(iset, ispin), Dthrs(iset, ispin)    
!             write (*,*) nctrs(iset, ispin), nspin
              if (istat.ne.0) then
                ifinputok=.false.
                goto 2233
              end if
!             write (54,*) nctrs(iset, ispin), ONthrs(iset, ispin), Dthrs(iset, ispin)
            end do
          end do


!         write (54,*) '----------------------------------------------------------'
 
          close (14)

2233    END SUBROUTINE read_input


        SUBROUTINE read_coords(vol, lvecs, atcoords, nat, ifdirect)

          IMPLICIT NONE
          integer, intent (in) :: nat
          real*8, intent (out) :: vol
          real*8, intent (out), dimension (3,3) :: lvecs (3,3)
          real*8, intent (out), dimension (inp%natom,3) :: atcoords
          logical, intent (out) :: ifdirect
          character :: coordtype
          integer :: iat, icoord, ivec, idmn

          ifdirect=.true.
          open (19, file='CONTCAR', status='old')
          read (19,*)
          read (19,*), vol
          do ivec=1,3
            read(19,*), (lvecs(ivec,idmn), idmn=1,3)
          end do
          read (19,*)
          read (19,*), coordtype

          if (coordtype=='c'.or.coordtype=='C'.or.coordtype=='k'.or.coordtype=='K') then
            ifdirect=.false.
          end if

          do iat=1, nat
            read(19,*), (atcoords(iat,icoord), icoord=1,3)
          end do
    
        END SUBROUTINE read_coords


      
   END SUBROUTINE SSAdNDP

END MODULE SSAdNDP_module
