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

! M_psi.f90
! Module Description
! ===========================================================================
!     This module is responsible for solving the Schrodinger equation. The
! main subroutine in this module to call is psirc, which will solve the
! Schrodinger equation for a given shell with a certain potential.
!
!       calculate_psi.f90 - this routine will find the solution of the Schrodinger
! equation for a given potential (vee + vxc + ...)
!       integrate_hpsi.f90 - integrates Schordinger equation - we are solving
! energy to then iterate to get new psi in get_psi.f90
!       get_psi.f90 - after integration - find now new psi for given energy
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
        module M_psi

! Type Declaration
! ===========================================================================

! Variable Declaration and Description
! ===========================================================================

! module procedures
        contains


! ===========================================================================
! calculate_psi
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine finds the wavefunction and energy for a potential given
! the boundary condition that psi = 0 at rcutoff. (Hartree, Bohr units). This
! routine can do any l state. The eigenvalue is negative for a bound state.
! This routine will solve for any excited state, use nexcite = 0 for the usual
! node_loweress case. Non-zero nexcit will give the corresponding excited
! state with angular momentum l.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine calculate_psi (meshpsi, nexcite, l, rcutoff, v, eout, psi)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: meshpsi
        integer, intent (in) :: l
        integer, intent (in) :: nexcite

        real, intent (in) :: rcutoff
        real, intent (in), dimension (meshpsi) :: v

! Output
        real, intent (out) :: eout
        real, intent (out), dimension (meshpsi) :: psi

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: small = 1.0d-14

! Local Variable Declaration and Description
! ===========================================================================
        integer iteration
        integer node
        integer node_lower, node_upper

        real alnd, alnd1, alnd2
        real alnd_lower, alnd_upper
        real ebind, ebind1, ebind2
        real ebind_lower, ebind_upper

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Assume all energies are not orders of magnitude away from atomic units --
! otherwise you might want to put a guessed energy here to start.
! Find upper bound to binding energy
        ebind_upper = 1.0d0
        do iteration = 1, 100
          ebind = 2.0d0*ebind_upper
          call integrate_hpsi (meshpsi, l, rcutoff, ebind, v, node, alnd)
          if ((node .lt. nexcite) .or. (node .eq. nexcite .and. alnd .gt. 0.0d0)) exit
          ebind_upper = ebind
        end do
        if (.not. ((node .lt. nexcite) .or. (node .eq. nexcite .and. alnd .gt. 0.0d0))) stop ' no upper bound found -- stopping '
        node_upper = node
        alnd_upper = alnd
        ebind_upper = ebind

! Find lower bound to binding energy
        ebind_lower = -1.0d0
        do iteration = 1, 100
          ebind = 2.0d0*ebind_lower
          call integrate_hpsi (meshpsi, l, rcutoff, ebind, v, node, alnd)
          if (.not. ((node .lt. nexcite) .or. (node .eq. nexcite .and. alnd .gt. 0.0d0))) exit
          ebind_lower = ebind
        end do
        if ((node .lt. nexcite) .or. (node .eq. nexcite .and. alnd .gt. 0.0d0)) stop ' no upper bound found -- stopping '
        node_lower = node
        alnd_lower = alnd
        ebind_lower = ebind

! Do binary chop to get close
        do iteration = 1, 100
          if (node_lower .eq. node_upper) exit
          ebind = 0.5d0*(ebind_upper + ebind_lower)
          call integrate_hpsi (meshpsi, l, rcutoff, ebind, v, node, alnd)
          if ((node .lt. nexcite) .or. (node .eq. nexcite .and. alnd .gt. 0.0d0)) then
            ebind_upper = ebind
            node_upper = node
            alnd_upper = alnd
          else
            ebind_lower = ebind
            node_lower = node
            alnd_lower = alnd
          end if
        end do
        if (node_lower .ne. node_upper) then
          write (*,*) ' ******************* WARNING ************************* '
          write (*,*) ' No convergence in binary chop. Generally, this means  '
          write (*,*) ' your grid is too course or something else is causing  '
          write (*,*) ' an extra node. '
        end if

! Linearly interpolate to get good eigenvalues
        ebind1 = ebind_upper
        ebind2 = ebind_lower
        alnd1 = alnd_upper
        alnd2 = alnd_lower
        do iteration = 1, 100
          if (abs(alnd1-alnd2) .lt. small) exit
          ebind = ebind1 - (ebind2 - ebind1)*alnd1/(alnd2 - alnd1)
          if (ebind .gt. ebind_upper .or. ebind .lt. ebind_lower)          &
     &      ebind = 0.5d0*(ebind_upper + ebind_lower)
          call integrate_hpsi (meshpsi, l, rcutoff, ebind, v, node, alnd)
          if ((node .lt. nexcite) .or. (node .eq. nexcite .and. alnd .gt. 0.0d0)) then
            ebind_upper = ebind
            alnd_upper = alnd
          else
            ebind_lower = ebind
            alnd_lower = alnd
          end if
          alnd2 = alnd1
          alnd1 = alnd
          ebind2 = ebind1
          ebind1 = ebind
        end do
        if (abs(alnd1-alnd2) .gt. small)  stop ' No convergence in linear extrapolation '
        eout = -ebind
        call get_psi (meshpsi, l, rcutoff, ebind, v, psi)

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
        return
        end subroutine calculate_psi


! ===========================================================================
! integrate_hpsi
! ===========================================================================
! Subroutine Description
! ===========================================================================
!    This routines integrates the Schroedinger equation for a trial binding
! energy. The number of nodes are returned and the discontinuity in the
! logarithmic derivative.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine integrate_hpsi (meshpsi, l, rcutoff, ebind, v, node, alnd)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: meshpsi
        integer, intent (in) :: l

        real, intent (in) :: ebind
        real, intent (in) :: rcutoff
        real, intent (in), dimension (meshpsi) :: v

! Output
        integer, intent (out) :: node

        real, intent (out) :: alnd

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer imatch
        integer ipoint
        integer istart

        real dr
        real h12
        double precision :: u1, u2, u3
        real v1, v2, v3

        real, dimension (:), allocatable :: r
        real, dimension (:), allocatable :: v_ang

! Allocate Arrays
! ===========================================================================
        allocate (r (meshpsi))
        allocate (v_ang (meshpsi))

! Procedure
! ===========================================================================
        dr = rcutoff/(meshpsi - 1)
        h12 = dr**2/12.0d0
        do ipoint = 2, meshpsi
          r(ipoint) = (ipoint - 1)*dr
        end do
        v_ang(1) = 0.0d0
        v_ang(2:meshpsi) = l*(l + 1.0d0)/r(2:meshpsi)**2

! imatch = matching point for the integration
        istart = 1
        imatch = meshpsi/2

! Initialize node counter
        node = 0

! Integrate out from origin
        u1 = 0.0d0
        u2 = dr
        v1 = 0.0d0
        v2 = v(istart+1) + v_ang(istart+1) + ebind
        do ipoint = istart + 2, imatch
          v3 = v(ipoint) + v_ang(ipoint) + ebind
          u3 = (u2*(2.0d0 + 10.0d0*h12*v2) - u1*(1.0d0 - h12*v1))/(1.0d0 - h12*v3)
          if (sign(1.0d0,u3) .ne. sign(1.0d0,u2)) node = node + 1
          u1 = u2
          u2 = u3
          v1 = v2
          v2 = v3
        end do
        v3 = v(imatch+1) + v_ang(imatch+1) + ebind
        u3 = (u2*(2.0d0 + 10.0d0*h12*v2) - u1*(1.0d0 - h12*v1))/(1.0d0 - h12*v3)
        alnd = (u3 - u1)/(2.0d0*dr*u2)

! Now integrate in from rmax
        u1 = 0.0d0
        u2 = dr
        v1 = v(meshpsi) + v_ang(meshpsi) + ebind
        v2 = v(meshpsi-1) + v_ang(meshpsi-1) + ebind
        do ipoint = meshpsi-2, imatch, -1
          v3 = v(ipoint) + v_ang(ipoint) + ebind
          u3 = (u2*(2.0d0 + 10.0d0*h12*v2) - u1*(1.0d0 - h12*v1))            &
     &         /(1.0d0 - h12*v3)
          if (sign(1.0d0,u3) .ne. sign(1.0d0,u2)) node = node + 1
          u1 = u2
         u2 = u3
         v1 = v2
         v2 = v3
       end do
       v3 = v(meshpsi-imatch-1) + v_ang(meshpsi-imatch-1) + ebind
       u3 = (u2*(2.0d0 + 10.0d0*h12*v2) - u1*(1.0d0 - h12*v1))/(1.0d0 - h12*v3)

! Calculate discontinuity in log derivative - it may look like this is lower
! order, but matching these actually matches the solution exactly
      alnd = alnd + (u3 - u1)/(2.0d0*dr*u2)

! Deallocate Arrays
! ===========================================================================
      deallocate (r)
      deallocate (v_ang)

! Format Statements
! ===========================================================================
      return
      end subroutine integrate_hpsi


! ===========================================================================
! get_psi
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routines integrates the Schroedinger equation for a trial
! binding energy. The number of nodes are returned and the discontinuity in
! the logarithmic derivative.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine get_psi (meshpsi, l, rcutoff, ebind, v, psi)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: meshpsi
        integer, intent (in) :: l

        real, intent (in) :: ebind
        real, intent (in) :: rcutoff
        real, intent (in), dimension (meshpsi) :: v

! Output
        real, intent (out), dimension (meshpsi) :: psi

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer imatch
        integer ipoint
        integer istart

        real anorm
        real dr
        real h12
        real ratio
        real u1, u2, u3
        real v1, v2, v3

        real, dimension (:), allocatable :: r
        real, dimension (:), allocatable :: v_ang

! Allocate Arrays
! ===========================================================================
        allocate (r (meshpsi))
        allocate (v_ang (meshpsi))

! Procedure
! ===========================================================================
        dr = rcutoff/(meshpsi - 1)
        h12 = dr**2/12.0d0
        do ipoint = 2, meshpsi
          r(ipoint) = (ipoint - 1)*dr
        end do
        v_ang(1) = 0.0d0
        v_ang(2:meshpsi) = l*(l + 1.0d0)/r(2:meshpsi)**2

! imatch = matching point for the integration
        istart = 1
        imatch = meshpsi/2

! Set up initial solution at origin and rmax
        psi(1) = 0.0d0
        psi(2) = dr
        psi(meshpsi) = 0.0d0
        psi(meshpsi-1) = dr

! Integrate out from origin
        u1 = 0.0d0
        u2 = dr
        v1 = v(istart) + ebind
        v2 = v(istart+1) + v_ang(istart+1) + ebind
        do ipoint = istart + 2, imatch
          v3 = v(ipoint) + v_ang(ipoint) + ebind
          u3 = (u2*(2.0d0 + 10.0d0*h12*v2) - u1*(1.d0 - h12*v1))/(1.0d0 - h12*v3)
          psi(ipoint) = u3
          u1 = u2
          u2 = u3
          v1 = v2
          v2 = v3
        end do
        ratio = u3

! Now integrate in from rmax
        u1 = 0.0d0
        u2 = dr
        v1 = v(meshpsi) + v_ang(meshpsi) + ebind
        v2 = v(meshpsi-1) + v_ang(meshpsi-1) + ebind
        do ipoint = meshpsi-2, imatch, -1
          v3 = v(ipoint) + v_ang(ipoint) + ebind
          u3 = (u2*(2.0d0 + 10.0d0*h12*v2) - u1*(1.0d0 - h12*v1))/(1.0d0 - h12*v3)
          psi(ipoint) = u3
          u1 = u2
          u2 = u3
          v1 = v2
          v2 = v3
        end do
        ratio = ratio/u3

! Match solutions
        psi(imatch:meshpsi) = psi(imatch:meshpsi)*ratio

! normalize using trapezoidal rule -- remember the end points are zero
        anorm = 1.0d0/sqrt(dot_product(psi,psi)*dr)
        psi = anorm*psi

! Deallocate Arrays
! ===========================================================================
        deallocate (r)
        deallocate (v_ang)

! Format Statements
! ===========================================================================
! None

        return
        end subroutine get_psi

! End Module
! ===========================================================================
        end module M_psi
