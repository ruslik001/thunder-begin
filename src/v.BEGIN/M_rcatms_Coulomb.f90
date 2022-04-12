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

! M_rcatms_Coulomb
! Module Description
! ===========================================================================
!		Module that creates the fireball wavefunctions - radial part.
! This module contains the following routines:
!
!       calculate_vee.f90 - calculates the Coulomb potential
!       calculate_vxc.90 - calculates the exchange-correlation potential
!       calculate_exc.f90 - calculates the exchange-correlation energy
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
        module M_rcatms_Coulomb
        use M_atom_functions
        use M_atomPP_functions
        use M_vxc_Harris

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains


! ===========================================================================
! calculate_vee
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This subroutine evaluates the electron-electron Hartree potential for
! a given density. It is given by:
!
!    vee(r) = 2*int(dr' n(r)/!r-r'!)
!           = 2*(1/r*int[0,r] sigma(t)dt + int[r,rc] sigma(t)/t dt  )
!
! Let x1(r) = int[0,r] sigma(t)dt
!     x2(r) = int[r,rc] sigma(t)/t dt.
!
! then, vee(r) = 2/r*x1(r) + 2*(x2(wfmesh) - x2(i))
!
! Integrate to get x1(r), via trapezoidal rule.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine calculate_vee (ispecies)
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
        integer mesh                     ! mesh size for the wf

        real a1, a2
        real b1, b2
        real dr                          ! distance between grid points
        real r                           ! value of radial point

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize some quantities:
        mesh = wf(ispecies)%mesh_max
        dr = wf(ispecies)%dr_min

! Find x1(r)
        a1 = 0.0d0
        wf(ispecies)%vee = 0.0d0
        r = -dr
        r = r + dr
        do ipoint = 2, mesh
          r = r + dr
          a2 = 0.5d0*wf(ispecies)%sigma(ipoint)
          a1 = a1 + a2
          wf(ispecies)%vee(ipoint) = 2.0d0*a1*dr/r
          a1 = a1 + a2
        end do

! Next get x2(r) - also trapezoidal rule
        b1 = 0.0d0
        do ipoint = mesh, 2, -1
          b2 = 0.5d0*wf(ispecies)%sigma(ipoint)/r
          b1 = b1 + b2
          wf(ispecies)%vee(ipoint) = wf(ispecies)%vee(ipoint) + 2.0d0*b1*dr
          b1 = b1 + b2
          r = r - dr
        end do
        wf(ispecies)%vee(1) = 2.0d0*b1*dr

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None
        return
        end subroutine calculate_vee


! ===========================================================================
! calculate_vxc
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This function computes the potential vxc with the densities ni of atom in_i.
!
! On input:  in1: atomic indices
!            r, z: geometry information for the charge gradient
!            ix: switch for the derivatives
!
! On output:  vxc: vxc[n1(r1)]  (the one center vxc)

! We calculate vxc(n1). We compute neutral atoms only here.
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
! Program Declaration
! ===========================================================================
        subroutine calculate_vxc (ispecies)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispecies  !< two centers species

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iexc                      !< which exchange-correlation
        integer ipoint                    !< loop over mesh points

        real dr                           !< distance between mesh points
        real r                            !< value of r in Bohr radii

! Value of density and corresponding derivatives at the point r, z
!        real rescale                      !< rescale densities
        real density
        real density_p, density_pp

! Output from calling get_potxc_1c
        real exc
        real vxc
        real dnuxc
        real dnuxcs
        real dexc

        real xc_fraction                  !< fraction of exact exchange

! Procedure
! ===========================================================================
! Initialize some variables
        iexc = species_PP(ispecies)%iexc
        xc_fraction = species_PP(ispecies)%xc_fraction

! Initialize dr for this one-center case.
        dr = wf(ispecies)%dr_min
        mesh = wf(ispecies)%mesh_max

! One-center piece: vxc[n1(r1)]
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Evaluate the density for the one-center - in1
        r = - dr
        r = r + dr
        do ipoint = 2, mesh - 1
          r = r + dr
          density = wf(ispecies)%rho(ipoint)
          density_p = (wf(ispecies)%rho(ipoint + 1)                          &
     &                 - wf(ispecies)%rho(ipoint - 1))/(2.0d0*dr)
          density_pp = (wf(ispecies)%rho(ipoint + 1)                         &
     &                  - 2.0d0*wf(ispecies)%rho(ipoint) + wf(ispecies)%rho(ipoint - 1))/dr**2
          call get_potxc_1c (iexc, xc_fraction, r, density, density_p,       &
     &                       density_pp, exc, vxc, dnuxc, dnuxcs, dexc)
          wf(ispecies)%vxc(ipoint) = 2.0d0*vxc
        end do

! endpoints of vxc
        ! ipoint = mesh
        density = wf(ispecies)%rho(mesh)
        density_p = (2.0d0*(wf(ispecies)%rho(mesh) - wf(ispecies)%rho(mesh-2)) &
     &               - (wf(ispecies)%rho(mesh-1) - wf(ispecies)%rho(mesh-3)))/(2.0d0*dr)
        density_pp = (2.0d0*(wf(ispecies)%rho(mesh)                          &
     &                       - 2.0d0*wf(ispecies)%rho(mesh-1)                &
     &                       + wf(ispecies)%rho(mesh-2))                     &
     &                - (wf(ispecies)%rho(mesh-1)                            &
     &                   - 2.0d0*wf(ispecies)%rho(mesh-2)                    &
     &                   + wf(ispecies)%rho(mesh-3)))/dr**2
        call get_potxc_1c (iexc, xc_fraction, r, density, density_p,         &
     &                     density_pp, exc, vxc, dnuxc, dnuxcs, dexc)
        wf(ispecies)%vxc(mesh) = 2.0d0*vxc

        ! ipoint = 1
        r = 0.0d0
        density = wf(ispecies)%rho(1)
        density_p = (2.0d0*(wf(ispecies)%rho(3) - wf(ispecies)%rho(1))       &
     &               - (wf(ispecies)%rho(4) - wf(ispecies)%rho(2)))/(2.0d0*dr)
        density_pp = (2.0d0*(wf(ispecies)%rho(3) - 2.0d0*wf(ispecies)%rho(2) &
     &                       + wf(ispecies)%rho(1))                          &
     &                - (wf(ispecies)%rho(4) - 2.0d0*wf(ispecies)%rho(3)     &
     &                   + wf(ispecies)%rho(2)))/dr**2
        call get_potxc_1c (iexc, xc_fraction, r, density, density_p,         &
     &                     density_pp, exc, vxc, dnuxc, dnuxcs, dexc)
        wf(ispecies)%vxc(1) = 2.0d0*vxc

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine calculate_vxc


! ===========================================================================
! calculate_exc
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This function computes only exc(n1) with the densities ni of atom in_i.
!
! On input:  in1: atomic indices
!            r, z: geometry information for the charge gradient
!            ix: switch for the derivatives
!
! On output:  exc: exc[n1(r1)]  (the one center exc)

! We calculate exc(n1). We compute neutral atoms only here.
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
! Program Declaration
! ===========================================================================
        subroutine calculate_exc (ispecies)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispecies  !< two centers species

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iexc                      !< which exchange-correlation
        integer ipoint                    !< loop over mesh points

        real dr                           !< distance between mesh points
        real r                            !< value of r in Bohr radii
!       real rmax                         !< maximum value of radial grid

! Value of density and corresponding derivatives at the point r, z
!       real rescale                      !< rescale densities
        real density
        real density_p, density_pp

! Output from calling get_potxc_1c
        real exc
        real vxc
        real dnuxc
        real dnuxcs
        real dexc

        real xc_fraction                  !< fraction of exact exchange

! Procedure
! ===========================================================================
! Initialize some variables
        iexc = species_PP(ispecies)%iexc
        xc_fraction = species_PP(ispecies)%xc_fraction

! Initialize dr for this one-center case.
        dr = wf(ispecies)%dr_min
        mesh = wf(ispecies)%mesh_max

! One-center piece: exc[n1(r1)]
! Compute the exchange correlation energy for the one-center case
! ***************************************************************************
        wf(ispecies)%exc = 0.0d0
        r = - dr
        r = r + dr
        do ipoint = 2, mesh - 1
          r = r + dr
          density = wf(ispecies)%rho(ipoint)
          density_p = (wf(ispecies)%rho(ipoint + 1)                          &
     &                 - wf(ispecies)%rho(ipoint - 1))/(2.0d0*dr)
          density_pp = (wf(ispecies)%rho(ipoint + 1)                         &
     &                  - 2.0d0*wf(ispecies)%rho(ipoint) + wf(ispecies)%rho(ipoint - 1))/dr**2
          call get_potxc_1c (iexc, xc_fraction, r, density, density_p,       &
     &                       density_pp, exc, vxc, dnuxc, dnuxcs, dexc)
          wf(ispecies)%exc =                                                 &
     &      wf(ispecies)%exc + dr*r**2*wf(ispecies)%rho(ipoint)*(exc - vxc)
        end do

! endpoints of vxc
        ! ipoint = mesh
        density = wf(ispecies)%rho(mesh)
        density_p = (2.0d0*(wf(ispecies)%rho(mesh) - wf(ispecies)%rho(mesh-2)) &
     &               - (wf(ispecies)%rho(mesh-1) - wf(ispecies)%rho(mesh-3)))/(2.0d0*dr)
        density_pp = (2.0d0*(wf(ispecies)%rho(mesh)                          &
     &                       - 2.0d0*wf(ispecies)%rho(mesh-1)                &
     &                       + wf(ispecies)%rho(mesh-2))                     &
     &                - (wf(ispecies)%rho(mesh-1)                            &
     &                   - 2.0d0*wf(ispecies)%rho(mesh-2)                    &
     &                   + wf(ispecies)%rho(mesh-3)))/dr**2
        call get_potxc_1c (iexc, xc_fraction, r, density, density_p,         &
     &                     density_pp, exc, vxc, dnuxc, dnuxcs, dexc)
        wf(ispecies)%exc =                                                   &
     &    wf(ispecies)%exc + dr*r**2*wf(ispecies)%rho(mesh)*(exc - vxc)

        ! ipoint = 1
        r = 0.0d0
        density = wf(ispecies)%rho(1)
        density_p = (2.0d0*(wf(ispecies)%rho(3) - wf(ispecies)%rho(1))       &
     &               - (wf(ispecies)%rho(4) - wf(ispecies)%rho(2)))/(2.0d0*dr)
        density_pp = (2.0d0*(wf(ispecies)%rho(3) - 2.0d0*wf(ispecies)%rho(2) &
     &                       + wf(ispecies)%rho(1))                          &
     &                - (wf(ispecies)%rho(4) - 2.0d0*wf(ispecies)%rho(3)     &
     &                   + wf(ispecies)%rho(2)))/dr**2
        call get_potxc_1c (iexc, xc_fraction, r, density, density_p,         &
     &                     density_pp, exc, vxc, dnuxc, dnuxcs, dexc)
        wf(ispecies)%exc =                                                   &
     &    wf(ispecies)%exc + dr*r**2*wf(ispecies)%rho(1)*(exc - vxc)

        ! rescale by 4*pi
        wf(ispecies)%exc = 2.0d0*(4.0d0*4.0d0*atan(1.0d0))*wf(ispecies)%exc

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine calculate_exc


! End Module
! ===========================================================================
        end module M_rcatms_Coulomb
