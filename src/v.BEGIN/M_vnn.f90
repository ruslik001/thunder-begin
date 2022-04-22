! copyright info:
!
!                             @Copyright 2022
!                           Fireball Committee
! Hong Kong Quantum AI Laboratory, Ltd. - James P. Lewis, Chair
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek
! Arizona State University - Otto F. Sankey

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! California Institute of Technology - Brandon Keith
! Czech Institute of Physics - Prokop Hapala
! Czech Institute of Physics - Vladimír Zobač
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
! Synfuels China Technology Co., Ltd. - Pengju Ren
! Washington University - Pete Fedders
! West Virginia University - Ning Ma and Hao Wang
! also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
!      and Spencer Shellman

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! M_vnn.f90
! Module Description
! ===========================================================================
!     This module is responsible for creating the Hartree potentials and
! write them out to files (called Z.na0, Z.na1, Z.ena1, etc.).
!
!       calculate_vnn.f90 - here the Hartree potential is calculated
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!

! Module Declaration
! ===========================================================================
        module M_vnn
        use M_species
        use M_atom_functions
        use M_atomPP_functions
        use M_atomPP_ion_functions

! Type Declaration
! ===========================================================================

! Variable Declaration and Description
! ===========================================================================

! module procedures
        contains


! ===========================================================================
! initialize_na.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine will initialize and populate some of the arrays for
! the type na.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine initialize_na
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: mesh_na_start = 801

! Local Variable Declaration and Description
! ===========================================================================
!       integer ipoint                    ! loop over mesh points
        integer ispecies                  ! loop over species
        integer issh                      ! loop over shells
        integer mesh                      ! new neutral atom potential mesh
        integer nssh                      ! loop over maximum shells - nssh

!       real r, dr                        ! value of grid point and dr
        real dr

! Allocate Arrays
! ===========================================================================
        allocate (na (nspecies))

! Procedure
! ===========================================================================
        do ispecies = 1, nspecies
          nssh = species(ispecies)%nssh

          allocate (na(ispecies)%shell_data(0:nssh))
          na(ispecies)%dr_min = species(ispecies)%rcutoffA_max/real(mesh_na_start - 1)
          na(ispecies)%mesh_max = mesh_na_start

! Allocate the arrays according to mesh.
          mesh = na(ispecies)%mesh_max
          dr =  na(ispecies)%dr_min

          allocate (na(ispecies)%r(mesh))
          do issh = 0, species(ispecies)%nssh
            allocate (na(ispecies)%shell_data(issh)%FofR(mesh))
          end do
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine initialize_na


! ===========================================================================
! calculate_vnn.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine calculates the vnn potential and shell potentials.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine calculate_vnn (ispecies)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer ispecies                 ! which species are we calculating

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint, jpoint            ! loop over mesh points
        integer issh                      ! loop over shells
!       integer iten, ione, itwo          ! character array places
        integer iten, ione                ! character array places
        integer lssh                      ! quantum number l for shell
        integer nssh                      ! loop over maximum shells - nssh
        integer nZ                        ! atomic number

        integer mesh                      ! new neutral atom potential mesh

        real alpha                        ! alpha of the pesudopotential
        real factor                       ! integration factor

        real pi

        real r1, r2, dr1, dr2             ! integration points and grid spacings

        real r, dr                        ! value of grid point and dr
        real rmax                         ! maximum value of r along the grid

        real rcutoff                      ! cutoff radius
        real vcore                        ! nuclear core potential

        real, dimension (:), allocatable :: term1
        real, dimension (:), allocatable :: term2

        real, dimension (:), allocatable :: xnocc    ! occupation number

        ! for writing out wavefunctions
!       integer inum, iremainder
        real, dimension (:), allocatable :: xx, yy

        character (len = 3) buffer
        character (len = 30) filename

        character (len = 1), dimension (0:9) :: z

        logical formatted

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize pi
        pi = 4.0d0*atan(1.0d0)

! Initialize some symbols
        do ipoint = 0, 9
          z(ipoint) = char(48 + ipoint)
        end do

        write (ilogfile,*)
        write (ilogfile,*) ' Beginning calculate_vnn subroutine for ispecies = ', ispecies
        write (ilogfile,*)

! Find our effective values for xnocc here based on rcatms_excite
        nssh = species(ispecies)%nssh
        allocate (xnocc (nssh))
        allocate (term1 (nssh))
        allocate (term2 (nssh))
        do issh = 1, nssh
          xnocc(issh) = species(ispecies)%shell(issh)%Qneutral
        end do

        dr =  na(ispecies)%dr_min
        mesh = na(ispecies)%mesh_max
        rcutoff = species(ispecies)%rcutoffA_max

        r = -dr
        do ipoint = 1, mesh
          r = r + dr
          term1 = 0.0d0
          term2 = 0.0d0
          na(ispecies)%r(ipoint) = r

          r1 = 0.0d0
          r2 = r
          dr1 = r/real(wf(ispecies)%mesh_max - 1)
          dr2 = (rcutoff - r)/real(wf(ispecies)%mesh_max - 1)
          do jpoint = 1, wf(ispecies)%mesh_max
            ! set up integration factors
            factor = 4.0d0
            if (mod(jpoint,2) .ne. 0) factor = 2.0d0
            if (jpoint .eq. 1 .or. jpoint .eq. wf(ispecies)%mesh_max) factor = 1.0d0
            do issh = 1, species(ispecies)%nssh
              rmax = species(ispecies)%shell(issh)%rcutoffA
              term1(issh) = term1(issh) + factor*(r1*psiofr(r1, rmax, ispecies, issh))**2
              term2(issh) = term2(issh) + factor*r2*psiofr(r2, rmax, ispecies, issh)**2
            end do
            r1 = r1 + dr1
            r2 = r2 + dr2
          end do

          do issh = 1, species(ispecies)%nssh
            term1(issh) = (dr1/3.0d0)*term1(issh)
            term2(issh) = (dr2/3.0d0)*term2(issh)

            na(ispecies)%shell_data(issh)%FofR(ipoint) = term2(issh)
            if (r .ne. 0.0d0) na(ispecies)%shell_data(issh)%FofR(ipoint) =   &
     &                        na(ispecies)%shell_data(issh)%FofR(ipoint) + term1(issh)/r
           end do

          alpha = species_PP(ispecies)%alpha
          if (r .ne. 0.0d0) then
            vcore = - (species(ispecies)%Zval/r)*erf(alpha*r)
          else
            vcore = - species(ispecies)%Zval*2.0d0*alpha/sqrt(pi)
          end if
          na(ispecies)%shell_data(0)%FofR(ipoint) = vcore + vPP_shortofr(r, ispecies)/P_eq2
          do issh = 1, species(ispecies)%nssh
            na(ispecies)%shell_data(0)%FofR(ipoint) =                        &
     &      na(ispecies)%shell_data(0)%FofR(ipoint)                          &
     &       + xnocc(issh)*na(ispecies)%shell_data(issh)%FofR(ipoint)
          end do

! Write out the information to the files
          iten = int(float(mesh)/10.0d0)
          if (mod(ipoint,iten) .eq. 1) then
            write (ilogfile,*)
            write (ilogfile,201) ipoint, r, na(ispecies)%shell_data(0)%FofR(ipoint)
            do issh = 1, species(ispecies)%nssh
              write (ilogfile,202) issh - 1,  na(ispecies)%shell_data(issh)%FofR(ipoint)
            end do
          end if
        end do

! ***************************************************************************
! Write Hartree potentials into output files
! ***************************************************************************
! ***************************************************************************
! Write wavefunction into output files
! ***************************************************************************
! First convert to angstrom units :    r   ---> r    /  (0.529177249)
!                                     r(r) ---> r(r) / (0.529177249)**1.5
! r(r) = u(r) / r
! Loop over the different states
        do issh = 0, species(ispecies)%nssh
          allocate (xx (mesh))
          allocate (yy (mesh))

          if (issh .eq. 0) then
            filename = species(ispecies)%na0file
          else
            filename = species(ispecies)%shell(issh)%nafile
          end if
          open (unit = 22, file = trim(Fdata_location)//'/'//trim(filename), &
     &          status = 'unknown', form = 'unformatted')
          if (issh .eq. 0) then
            write (22) species(ispecies)%rcutoff_max
          else
            write (22) species(ispecies)%shell(issh)%rcutoff,                &
     &                 species(ispecies)%rcutoff_max
          end if
         write (22) mesh
         if (issh .eq. 0) write (22) species(ispecies)%atomicE

! Find (approximately) the value at r=0
          xx(2) = na(ispecies)%r(2)
          xx(3) = na(ispecies)%r(3)

          yy(2) = na(ispecies)%shell_data(issh)%FofR(2)
          yy(3) = na(ispecies)%shell_data(issh)%FofR(3)

          xx(1) = 0.0d0
          yy(1) = - (yy(3) - yy(2))/(xx(3) - xx(2))*xx(2) + yy(2)

          xx(4:mesh) = na(ispecies)%r(4:mesh)
          yy(4:mesh) = na(ispecies)%shell_data(issh)%FofR(4:mesh)

          do ipoint = 1, mesh
            write (22) xx(ipoint), yy(ipoint)
          end do
          deallocate (xx)
          deallocate (yy)
          close (unit = 22)
        end do

! The variable symbol is a three character string which has the atomic number.
! For example, 001 is hydrogen, 008 is oxygen, and 094 is plutonium.
        nZ = species(ispecies)%nZ
        iten = nZ/10
        if (nZ .lt. 10) buffer(1:3) = z(0)//z(0)//z(nZ)
        if (nZ .ge. 10 .and. nZ .le. 99) then
          iten = nZ/10
          ione = nZ - 10*iten
          buffer(1:3) = z(0)//z(iten)//z(ione)
        end if

        do issh = 0, species(ispecies)%nssh
          if (issh .eq. 0) then
            write (filename,'(a3,".na.dat")') buffer(1:3)
            open (unit = 22, file = filename, status = 'unknown')
          else
            lssh = species(ispecies)%shell(issh)%lssh
            if (scan(species(ispecies)%shell(issh)%wffile,'e').ne.0) then
              if (lssh .eq. 0) then
                write (filename,'(a3,".na-s1.dat")') buffer(1:3)
                open (unit = 22, file = filename, status = 'unknown')
              else if (lssh .eq. 1) then
                write (filename,'(a3,".na-p1.dat")') buffer(1:3)
                open (unit = 22, file = filename, status = 'unknown')
              else if (lssh .eq. 2) then
                write (filename,'(a3,".na-d1.dat")') buffer(1:3) 
                open (unit = 22, file = filename, status = 'unknown')
              end if
            else
              if (lssh .eq. 0) then
                write (filename,'(a3,".na-s0.dat")') buffer(1:3)
                open (unit = 22, file = filename, status = 'unknown')
              else if (lssh .eq. 1) then
                write (filename,'(a3,".na-p0.dat")') buffer(1:3)
                open (unit = 22, file = filename, status = 'unknown')
              else if (lssh .eq. 2) then
                write (filename,'(a3,".na-d0.dat")') buffer(1:3)
                open (unit = 22, file = filename, status = 'unknown')
              end if
            end if
          end if
          do ipoint = 1, mesh
            write (22,101) na(ispecies)%r(ipoint), na(ispecies)%shell_data(issh)%FofR(ipoint)
          end do
          close (unit = 22)
        end do

        formatted = .false.
        if (formatted) then
! ***************************************************************************
! Write potential into FORMATTED output files
! ***************************************************************************
! Write out our vnn
          do issh = 0, species(ispecies)%nssh
            allocate (xx (mesh))
            allocate (yy (mesh))
    
            if (issh .eq. 0) then
              write (filename,'("formatted/",a14)') species(ispecies)%na0file(7:)
            else
              write (filename,'("formatted/",a14)') species(ispecies)%shell(issh)%nafile(7:)
            end if
            
            open (unit = 22, file = trim(Fdata_location)//'/'//trim(filename), &
  &               status = 'unknown')

            write (22, 301) filename(11:)
            write (22, *) species(ispecies)%nZ
            if (issh .eq. 0) then
              write (22, *) species(ispecies)%rcutoff_max
            else
              write (22, *) species(ispecies)%shell(issh)%rcutoff,             &
  &                         species(ispecies)%rcutoff_max
            end if
            write (22,*) mesh
            if (issh .eq. 0) write (22,*) species(ispecies)%atomicE

! Find (approximately) the value at r=0
            xx(2) = na(ispecies)%r(2)
            xx(3) = na(ispecies)%r(3)

            yy(2) = na(ispecies)%shell_data(issh)%FofR(2)
            yy(3) = na(ispecies)%shell_data(issh)%FofR(3)

            xx(1) = 0.0d0
            yy(1) = - (yy(3) - yy(2))/(xx(3) - xx(2))*xx(2) + yy(2)

            xx(4:mesh) = na(ispecies)%r(4:mesh)
            yy(4:mesh) = na(ispecies)%shell_data(issh)%FofR(4:mesh)

            do ipoint = 1, mesh
              write (22, 302) xx(ipoint), yy(ipoint)
            end do
            deallocate (xx)
            deallocate (yy)
            close (unit = 22)
          end do
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (2x, f6.3, 2x, f12.6)
201     format (2x, ' ipoint = ', i4, 2x, ' r = ', f9.5, 2x, ' V_na = ', d16.5)
202     format (22x, ' V_rho(issh = ', i1, ') =', d16.5)
301     format (2x, a11)
302     format (2d24.16)

        return
        end subroutine calculate_vnn


! ===========================================================================
! destroy_na
! ===========================================================================
! Subroutine Description
! ===========================================================================
!
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine destroy_na (nspecies)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nspecies

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ispecies

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! None

! Deallocate Arrays
! ===========================================================================
        do ispecies = 1, nspecies
          deallocate (na(ispecies)%r)
        end do

! Format Statements
! ===========================================================================
! None

        return
        end subroutine destroy_na


! End Module
! ===========================================================================



        end module M_vnn
