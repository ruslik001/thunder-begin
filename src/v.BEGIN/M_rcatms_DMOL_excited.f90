! copyright info:
!
!                             @Copyright 2013
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! Caltech - Brandon Keith
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
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

! M_rcatms_excited
! Module Description
! ===========================================================================
!		Module that creates the fireball wavefunctions - radial part.
! This module contains the following routines:
!
!       intializa_wf_ion.f90 - initialize the wavefunctions
!       calculate_rcatm_excited.f90 - calculate the wavefunction for cutoff
!       calculate_vconfine_excited.f90 - calculate the confinement potentials
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
        module M_rcatms_excited
        use M_species
        use M_atom_functions
        use M_atomPP_functions
        use M_atomPP_ion_functions
        use M_psi
        use M_rcatms_Coulomb_ion

! Type Declaration
! ===========================================================================
        real, allocatable :: vconfine_excited (:, :)   ! confinement potential

! module procedures
        contains


! ===========================================================================
! initialize_wf_ion.f90
! ===========================================================================
! Program Description
! ===========================================================================
!       Sets up the wfmesh grid for BEGIN by calculating the change in radius
! between wfmesh points and the radius at each wfmesh point.
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
! Declaration
! ===========================================================================
        subroutine initialize_wf_ion
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: mesh_wf_start = 1501

! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint                    ! loop over mesh points
        integer ispecies                  ! loop over species
        integer issh                      ! loop over shells
        integer mesh                      ! new wavefunction mesh
        integer lssh                      ! l quantum number for shell issh
        integer nssh                      ! number of shells

        real a0                           ! constant for radial wavefunction
        real r, dr                        ! value of grid point and dr
        real rcutoff                      ! cutoff radius
        real xnocc                        ! occupation number
        real xnorm                        ! normalization integral (= 1.0d0)

        real, dimension (:), allocatable :: integrand

        interface
          function simpson (mesh, integrand, dr)
            integer, intent (in) :: mesh
            real, intent (in), dimension (mesh) :: integrand
            real, intent (in) :: dr
            real simpson
          end function simpson
        end interface

! Allocate Arrays
! ===========================================================================
        allocate (wf_ion (nspecies))
        allocate (species_PP_ion (nspecies))

! Procedure
! ===========================================================================
! Loop over species and initialize the wavefunctions
        do ispecies = 1, nspecies
          nssh = species(ispecies)%nssh

          allocate (wf_ion(ispecies)%shell_data(nssh))
          wf_ion(ispecies)%dr_min = species(ispecies)%rcutoff_max/real(mesh_wf_start - 1)
          wf_ion(ispecies)%mesh_max = -1001

          ! The rcutoffs are different for each shell. Now re-evaluate what the
          ! true mesh should be for each shell.
          do issh = 1, species(ispecies)%nssh
            lssh = species(ispecies)%shell(issh)%lssh
            xnocc = species(ispecies)%shell(issh)%Qneutral_ion
            a0 = species(ispecies)%shell(issh)%a0

            ! find initial values of mesh points and dr
            dr = species(ispecies)%rcutoff_max/real(mesh_wf_start - 1)
            mesh = nint((species(ispecies)%shell(issh)%rcutoff)/dr) + 1
            wf_ion(ispecies)%shell_data(issh)%mesh = mesh
            wf_ion(ispecies)%mesh_max = max(wf_ion(ispecies)%mesh_max,       &
     &                                  wf_ion(ispecies)%shell_data(issh)%mesh)

            ! adjust slightly to make odd number of points and scale according
            ! to maximum cutoff
            rcutoff = real(mesh - 1)*dr
            species(ispecies)%shell(issh)%rcutoff = rcutoff
            dr = rcutoff/real(mesh - 1)
            wf_ion(ispecies)%shell_data(issh)%dr = dr

! Allocate the arrays according to new mesh.
            allocate (integrand (mesh))
            allocate (wf_ion(ispecies)%shell_data(issh)%r(mesh))
            allocate (wf_ion(ispecies)%shell_data(issh)%FofR(mesh))

            integrand = 0.0d0
            r = - dr
            do ipoint = 1, wf_ion(ispecies)%shell_data(issh)%mesh
              r = r + dr
              wf_ion(ispecies)%shell_data(issh)%r(ipoint) = r
              if (r .lt. rcutoff) integrand(ipoint) = r**(lssh + 1)*exp(-r/a0)
            end do

! Set wavefunction value (normalized)
            wf_ion(ispecies)%shell_data(issh)%FofR(:) =                      &
     &        integrand(:)/sqrt(dot_product(integrand,integrand)*dr)

! Check normalization
            write (ilogfile,*)
            write (ilogfile,*) ' Checking normalization [NORM(l) will not be 1, yet]'
            r = - dr
            do ipoint = 1, wf_ion(ispecies)%shell_data(issh)%mesh
              r = r + dr
              integrand(ipoint) = wf_ion(ispecies)%shell_data(issh)%r(ipoint)**2 &
      &                          *wf_ion(ispecies)%shell_data(issh)%FofR(ipoint)**2
            end do

            xnorm = simpson (mesh, integrand, dr)
            write (ilogfile,101) issh, xnorm
            deallocate (integrand)
          end do
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (2x, ' NORM (shell = ', i1, ') = ', f16.12)

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_wf_ion


! ===========================================================================
! calculate_rcatm_excited
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This will set up the potentials and self consistency for
! solving the Schrodinger equation and determine the wavefunction for a
! given rcutoff value.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine calculate_rcatm_excited (ispecies)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer ispecies                 ! which species are we calculating

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: max_iterations = 120

        real, parameter :: beta = 0.08d0
        real, parameter :: tolerance = 1.0d-9

! Local Variable Declaration and Description
! ===========================================================================
        integer iexc                     ! exchange-correlation option
        integer ipoint                   ! loop over mesh points
        integer issh                     ! loop over shells
!        integer iten, ione, itwo         ! character array places
        integer iten, ione               ! character array places
        integer inum, iremainder
        integer iteration                ! iterate to scf solution
        integer lssh                     ! quantum number l of the shell
        integer nssh                     ! loop over maximum shells - nssh
        integer nZ                         ! atomic number

        real difference                  ! integral sigma - Qneutral
        real pi

        real enew                        ! new band-structure energy
        real eold                        ! old band-structure energy
        real etot                        ! total energy

        real r, dr                       ! value of grid point and dr

        real vee                         ! temporary vee value
        real vxc                         ! temporary vxc value
        real xc_fraction                 ! mixing factor for exact exchange

        real, dimension (:), allocatable :: sigma_old
        real, dimension (:), allocatable :: xnocc

        ! total potential on the grid
        real, dimension (:), allocatable :: v

        ! for orthogonalization
        real, dimension (:), allocatable :: psi
        real, dimension (:), allocatable :: psi_DMOL

        ! for writing out wavefunctions
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
        write (ilogfile,*) ' Beginning calculate_rcatm_excited subroutine for ispecies = ', ispecies
        write (ilogfile,*)

! Find our effective values for xnocc here based on rcatms_excite
        nssh = species(ispecies)%nssh
        allocate (xnocc (nssh))
        do issh = 1, nssh
          xnocc(issh) = species(ispecies)%shell(issh)%Qneutral_ion
        end do

! Get the short-range and non-local psuedopotential contributions.
        iexc = species_PP(ispecies)%iexc
        xc_fraction = species_PP(ispecies)%xc_fraction
        call calculate_vPP_ion (ispecies)

! Calculate the confinement potential on the grid.
        call calculate_vconfine_excited (ispecies)

! Calculate the initial charge density
        dr = wf(ispecies)%dr_min
        mesh = wf(ispecies)%mesh_max
        allocate (wf_ion(ispecies)%r(mesh))
        allocate (wf_ion(ispecies)%rho(mesh))
        allocate (wf_ion(ispecies)%sigma(mesh))
        allocate (sigma_old(mesh))

        r = - dr
        wf_ion(ispecies)%sigma = 0.0d0
        do ipoint = 1, mesh
          r = r + dr
          wf_ion(ispecies)%r(ipoint) = r
          do issh = 1, nssh
            wf_ion(ispecies)%sigma(ipoint) = wf_ion(ispecies)%sigma(ipoint)  &
     &        + xnocc(issh)*psiofr_ion(r, ispecies, issh)**2
          end do
        end do

        ! Make sure that the sum of the charge densities is equal to the
        ! atomic number with a defined level of tolerance
        difference = abs(dr*sum(wf_ion(ispecies)%sigma) - real(sum(xnocc(1:nssh))))
        if (difference .gt. 1.0d-6) then
          write (ilogfile,*) ' ispecies = ', ispecies
          write (ilogfile,*) ' sum = ', dr*sum(wf_ion(ispecies)%sigma),      &
     &                       ' sum = ', real(sum(xnocc(1:nssh)))
          stop ' Check: Sum of charge density - not equal to Qneutral '
        end if

! Allocate some arrays for the potentials
        allocate (wf_ion(ispecies)%vee(mesh))
        allocate (wf_ion(ispecies)%vxc(mesh))

! Set up the density for exchange-correlation potentials
        wf_ion(ispecies)%rho(2:mesh - 1) = wf_ion(ispecies)%sigma(2:mesh - 1)&
     &                                    /wf_ion(ispecies)%r(2:mesh - 1)**2
        wf_ion(ispecies)%rho(1) = 2.0d0*wf_ion(ispecies)%rho(2) - wf_ion(ispecies)%rho(3)
        wf_ion(ispecies)%rho(mesh) = 2.0d0*wf_ion(ispecies)%rho(mesh - 1)    &
     &    - wf_ion(ispecies)%rho(mesh - 2)

! Calculate the electron-electron Hartree and the exchange-correlation       &
! potential on the grid.
        call calculate_vee_ion (ispecies)
        call calculate_vxc_ion (ispecies)

! Set the initial exact exchange potential to be 0
!       vex = 0.0d0

! Now solve the Schroedinger equation iteratively for each state in turn.
! Also, calculate core and valence charge density.
! ***************************************************************************
! ** Begin the self-consistency loop here **
! ***************************************************************************
! Allocate the potential to the mesh size of the maximum mesh.
        allocate (v(mesh))

        write (ilogfile, 101)
        write (ilogfile, 102)

        eold = -99.0d0    ! initial value for eold (the old energy)
        do iteration = 1, max_iterations
          r = - dr
          r = r + dr

! Store the current value of sigma and reset it
          sigma_old = wf_ion(ispecies)%sigma
          wf_ion(ispecies)%sigma = 0.0d0

! Begin a loop over all of the shells
          do issh = 1, species(ispecies)%nssh
            ! Coulomb potential
            vee = wf_ion(ispecies)%vee(1)
            ! exchange-correlation potential
            vxc = wf_ion(ispecies)%vxc(1)

            ! first point on the grid
            if (iexc .eq. 12) then
            else
              v(1) = 2.0d0*(species_PP_ion(ispecies)%vc(1)                   &
     &                      + species_PP_ion(ispecies)%vnl(issh,1)) + vee    &
     &                      + vxc + vconfine_excited(issh, 1)
            end if

! Remaining points on the grid:
            do ipoint = 2, mesh
              r = r + dr
              ! Coulomb potential
              vee = wf_ion(ispecies)%vee(ipoint)
              ! exchange-correlation potential
              vxc = wf_ion(ispecies)%vxc(ipoint)

! Sum up all the potential terms:
              if (iexc .eq. 12) then
              else
                v(ipoint) = 2.0d0*(species_PP_ion(ispecies)%vc(ipoint)       &
     &                             + species_PP_ion(ispecies)%vnl(issh,ipoint))&
     &                     + vee + vxc + vconfine_excited(issh,ipoint)
              end if
            end do

! Solve Schroedinger equation for potential v
            call calculate_psi (wf_ion(ispecies)%shell_data(issh)%mesh, 0,   &
     &                          species(ispecies)%shell(issh)%lssh,          &
     &                          species(ispecies)%shell(issh)%rcutoff, v,    &
     &                          wf_ion(ispecies)%shell_data(issh)%eigenvalue,&
     &                          wf_ion(ispecies)%shell_data(issh)%FofR)

! Add to charge denisty
            do ipoint = 1, wf_ion(ispecies)%shell_data(issh)%mesh
              wf_ion(ispecies)%sigma(ipoint) = wf_ion(ispecies)%sigma(ipoint)&
     &          + xnocc(issh)*wf_ion(ispecies)%shell_data(issh)%FofR(ipoint)**2
            end do
          end do   ! end loop over sh ells

! Mix in old and new sigma's
          wf_ion(ispecies)%sigma = (1.0d0 - beta)*wf_ion(ispecies)%sigma + beta*sigma_old

! Set up the density for exchange-correlation potentials
          wf_ion(ispecies)%rho(2:mesh - 1) = wf_ion(ispecies)%sigma(2:mesh - 1) &
     &       /(4.0d0*(4.0d0*atan(1.0d0))*wf_ion(ispecies)%r(2:mesh - 1)**2)
          wf_ion(ispecies)%rho(1) = 2.0d0*wf_ion(ispecies)%rho(2) - wf_ion(ispecies)%rho(3)
          wf_ion(ispecies)%rho(mesh) = 2.0d0*wf_ion(ispecies)%rho(mesh - 1)     &
     &       - wf_ion(ispecies)%rho(mesh - 2)

! Calculate the electron-electron Hartree and the exchange-correlation
! potential on the grid.
          call calculate_vee_ion (ispecies)
          call calculate_vxc_ion (ispecies)

! Evaluate the double counted electron-electron repulsion and the exchange-
! correlation energy.
! Do the integral using a trapezoidal rule.
          wf_ion(ispecies)%eee = 0.0d0
          do ipoint = 2, mesh - 1
            wf_ion(ispecies)%eee =                                             &
     &        wf_ion(ispecies)%eee + 0.5d0*dr*wf_ion(ispecies)%sigma(ipoint)*wf_ion(ispecies)%vee(ipoint)
          end do
          call calculate_exc_ion (ispecies)

! Calculate the exact exchange potential - FIXME
          if (iexc .eq. 12) then
            do issh = 1, species(ispecies)%nssh
!             call exchange ()
            end do
          end if

! Calculate the total energy:  Etot = sum(eigenvalue) - 1/2*Vee  -1/4*Vxc
! (Minus sign to convert binding energies into real energies)
! Sum over the eigenvalues:
          wf_ion(ispecies)%ebs = 0.0d0
          do issh = 1, species(ispecies)%nssh
           wf_ion(ispecies)%ebs =                                            &
     &       wf_ion(ispecies)%ebs + xnocc(issh)*wf_ion(ispecies)%shell_data(issh)%eigenvalue
          end do

! Find the total energy and the new energy
          etot = wf_ion(ispecies)%ebs - wf_ion(ispecies)%eee + wf_ion(ispecies)%exc
          enew = abs(wf_ion(ispecies)%ebs + etot)

! Check for self-consistency, then we can exit the loop
          difference = abs((enew - eold)/enew)
          write (ilogfile,103) iteration, etot, wf_ion(ispecies)%ebs,        &
     &      wf_ion(ispecies)%eee, wf_ion(ispecies)%exc, difference
          if (difference .le. tolerance) exit

! Set the old energy to the new energy, repeat
          eold = enew
        end do  ! End iteration loop

! ***************************************************************************
! ** End the self-consistency loop here **
! ***************************************************************************
        write (ilogfile,*)
        do issh = 1, species(ispecies)%nssh
          write (ilogfile,104) species(ispecies)%shell(issh)%lssh,           &
     &      xnocc(issh), wf_ion(ispecies)%shell_data(issh)%eigenvalue
        end do

        if (difference .gt. tolerance) then
          write (ilogfile,*)
          write (ilogfile,*) ' No self consistency after iterations = ', max_iterations
          stop
        end if

! Now orthogonalize the excited state wavefunctions with the ground state wavefunctions
        do issh = 1, species(ispecies)%nssh
           mesh = wf(ispecies)%shell_data(issh)%mesh
           allocate (psi (mesh))
           allocate (psi_DMOL (mesh))
           psi = wf(ispecies)%shell_data(issh)%FofR
           psi_DMOL = wf_ion(ispecies)%shell_data(issh)%FofR
           dr = wf(ispecies)%shell_data(issh)%dr
           call orthogonalize_psi (ilogfile, issh, mesh, psi, psi_DMOL, dr)
           wf_ion(ispecies)%shell_data(issh)%FofR = psi_DMOL
           deallocate (psi, psi_DMOL)
         end do

! ***************************************************************************
! Write wavefunction into output files
! ***************************************************************************
! First convert to angstrom units :    r   ---> r    /  (0.529177249)
!                                     r(r) ---> r(r) / (0.529177249)**1.5
! r(r) = u(r) / r
! Loop over the different states
        do issh = 1, nssh
          mesh = wf(ispecies)%shell_data(issh)%mesh
          allocate (xx (mesh))
          allocate (yy (mesh))

          ! Here we set the Qneutral_ion = 0.0d0 for writing out
          species(ispecies)%shell(issh)%Qneutral_ion = 0.0d0

          write (filename,'("basis/",a12)') species(ispecies)%shell(issh)%wffile_excited
          open (unit = 22, file = trim(Fdata_location)//'/'//trim(filename), &
     &          status = 'unknown', form = 'unformatted')
          write (22) species(ispecies)%shell(issh)%wffile_excited
          write (22) species(ispecies)%name, species(ispecies)%nZ
          write (22) species(ispecies)%shell(issh)%lssh,                     &
     &                   species(ispecies)%shell(issh)%Qneutral_ion
          write (22) species(ispecies)%shell(issh)%rcutoff,                  &
     &                   species(ispecies)%rcutoff_max
          write (22) mesh

! Find (approximately) the value at r=0
          xx(2) = wf(ispecies)%r(2)*P_abohr
          xx(3) = wf(ispecies)%r(3)*P_abohr

          yy(2) = wf_ion(ispecies)%shell_data(issh)%FofR(2)                  &
     &           /((P_abohr)**(1.50d0)*wf(ispecies)%r(2))
          yy(3) = wf_ion(ispecies)%shell_data(issh)%FofR(3)                  &
     &           /((P_abohr)**(1.50d0)*wf(ispecies)%r(3))

          xx(1) = 0.0d0
          yy(1) = - (yy(3) - yy(2))/(xx(3) - xx(2))*xx(2) + yy(2)

          xx(4:mesh-1) = wf(ispecies)%r(4:mesh-1)*P_abohr
          yy(4:mesh-1) = wf_ion(ispecies)%shell_data(issh)%FofR(4:mesh-1)    &
     &                  /((P_abohr)**(1.50d0)*wf(ispecies)%r(4:mesh-1))

          xx(mesh) = wf(ispecies)%r(mesh)*P_abohr
          yy(mesh) = 0.0d0

          do ipoint = 1, mesh
            write (22) xx(ipoint), yy(ipoint)
          end do
          close (unit = 22)

! Write out the wavefunctions for plotting purposes.
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

          lssh = species(ispecies)%shell(issh)%lssh
          if (lssh .eq. 0) then
            write (filename,'(a3,".wf-s1.dat")') buffer(1:3)
            open (unit = 22, file = filename, status = 'unknown')
          else if (lssh .eq. 1) then
            write (filename,'(a3,".wf-p1.dat")') buffer(1:3)
            open (unit = 22, file = filename, status = 'unknown')
          else if (lssh .eq. 2) then
            write (filename,'(a3,".wf-d1.dat")') buffer(1:3)
            open (unit = 22, file = filename, status = 'unknown')
          end if

          do ipoint = 1, mesh
            write (22, 301) wf(ispecies)%r(ipoint), yy(ipoint)
          end do
          close (unit = 22)

          deallocate (xx)
          deallocate (yy)
        end do


        formatted = .false.
        if (formatted) then
! ***************************************************************************
! Write wavefunction into FORMATTED output files
! ***************************************************************************
! First convert to angstrom units :    r   ---> r    /  (0.529177249)
!                                     r(r) ---> r(r) /
!                                     (0.529177249)**1.5
! r(r) = u(r) / r
! Loop over the different states
          do issh = 1, species(ispecies)%nssh
            mesh = wf(ispecies)%shell_data(issh)%mesh
            allocate (xx (mesh))
            allocate (yy (mesh))

            write (filename,'("formatted/",a15)')species(ispecies)%shell(issh)%wffile_excited
            open (unit = 22, file =trim(Fdata_location)//'/'//trim(filename), &
     &            status = 'unknown')
            write (22, 201) species(ispecies)%shell(issh)%wffile_excited
            write (22, 202) species(ispecies)%nZ, species(ispecies)%name
            write (22, 203) mesh
            write (22, 204) species(ispecies)%shell(issh)%rcutoff,&
     &        species(ispecies)%rcutoff_max
            write (22, 205) species(ispecies)%shell(issh)%lssh

! Find (approximately) the value at r=0
            xx(2) = wf(ispecies)%r(2)*P_abohr
            xx(3) = wf(ispecies)%r(3)*P_abohr

            yy(2) = wf_ion(ispecies)%shell_data(issh)%FofR(2)                    &
     &             /((P_abohr)**(1.50d0)*wf(ispecies)%r(2))
            yy(3) = wf_ion(ispecies)%shell_data(issh)%FofR(3)                    &
     &             /((P_abohr)**(1.50d0)*wf(ispecies)%r(3))

            xx(1) = 0.0d0
            yy(1) = - (yy(3) - yy(2))/(xx(3) - xx(2))*xx(2) + yy(2)

            xx(4:mesh-1) = wf(ispecies)%r(4:mesh-1)*P_abohr
            yy(4:mesh-1) = wf_ion(ispecies)%shell_data(issh)%FofR(4:mesh-1)      &
     &                    /((P_abohr)**(1.50d0)*wf(ispecies)%r(4:mesh-1))

            xx(mesh) = wf(ispecies)%r(mesh)*P_abohr
            yy(mesh) = 0.0d0

            inum = idint(real(mesh)/4.0d0)
            iremainder = mesh - (inum*4)

            do ipoint = 1, mesh - iremainder*4
              write (22, 400) yy(ipoint), yy(ipoint+1), yy(ipoint+2), yy(ipoint+3)
            end do

            if (iremainder .eq. 1) then
              write (22, 400) yy(mesh)
            else if (iremainder .eq. 2) then
              write (22, 400) yy(mesh - 1), yy(mesh)
            else if (iremainder .eq. 3) then
              write (22, 400) yy(mesh - 2), yy(mesh - 1), yy(mesh)
            end if
            close (unit = 22)

          deallocate (xx)
          deallocate (yy)
        end do
      end if

! Deallocate Arrays
! ===========================================================================
        deallocate (sigma_old)
        deallocate (v)
        deallocate (xnocc)

! Format Statements
! ===========================================================================
101     format (2x, ' iteration ', 2x, ' e(tot) ', 2x, ' eigenvalue ', 3x,   &
     &          ' u(e-e) ', 6x, ' u(x-c) ', 4x, ' difference ')
102     format (2x, 75('='))
103     format (5x, i3, 5x, 1e11.4, 4(2x,e12.5))
104     format (2x, ' l = ', i1, ' contains ', f9.6, ' electrons --- ' ,     &
     &          ' eigenvalue = ', f12.7)

201     format (2x, a12)
202     format (2x, i5, 2x, a11)
203     format (2x, i8)
204     format (2x, 2(2x,f7.4),2x,"0.0000")
205     format (2x, i3)

301     format (2x, f6.3, 2x, f12.6)

400     format (4d18.10)

! End Subroutine
! ===========================================================================
        return
        end subroutine calculate_rcatm_excited


! ===========================================================================
! calculate_vconfine_excited
! ===========================================================================
! Subroutine Description
! ===========================================================================
!     Calculates the x-confinement potential to make the wavefunction
! smoother near the cutoff point.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine calculate_vconfine_excited (ispecies)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispecies

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint                   ! loop over grid points
        integer issh                     ! loop over shells
        integer mesh                     ! mesh size for the wf

        real dr                          ! distance between grid points
        real r                           ! value of radial point
        real r0, v0                      ! confinement potential parameters
        real rcutoff                     ! radial cutoff

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize some quantities:
        mesh = wf(ispecies)%mesh_max
        dr = wf(ispecies)%dr_min

! Allocate the size of the confinement potential array
        allocate (vconfine_excited (species(ispecies)%nssh, mesh))
        vconfine_excited = 0.0d0
        if (species(ispecies)%ioptimize .eq. 1) then
          do issh = 1, species(ispecies)%nssh
            r0 = species(ispecies)%shell(issh)%r0_excited
            V0 = species(ispecies)%shell(issh)%V0_excited
            rcutoff = species(ispecies)%shell(issh)%rcutoff
            if (V0 .gt. 1.0d-4) then
              r = -dr
              do ipoint = 1, mesh
                r = r + dr
                if (r .gt. r0) then
                  vconfine_excited(issh, ipoint) =                           &
     &               V0*exp(-(rcutoff - r0)/(r - r0))/(rcutoff - r + 0.001d0)
                end if
              end do
            end if
          end do
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine calculate_vconfine_excited


! ===========================================================================
! destroy_rcatm_excited
! ===========================================================================
! Subroutine Description
! ===========================================================================
!
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine destroy_rcatm_excited (ispecies)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispecies

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
! None

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! None

! Deallocate Arrays
! ===========================================================================
        write(*,*) ispecies
        deallocate (wf_ion(ispecies)%r)
        deallocate (wf_ion(ispecies)%rho)
        deallocate (wf_ion(ispecies)%sigma)
        deallocate (wf_ion(ispecies)%vee)
        deallocate (wf_ion(ispecies)%vxc)

        deallocate (vconfine_excited)

! Format Statements
! ===========================================================================
! None

        return
        end subroutine destroy_rcatm_excited

! End Module
! ===========================================================================
        end module M_rcatms_excited
