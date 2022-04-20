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

! orthogonalize.f90
! Program Description
! ===========================================================================
!      This routine will orthogonalioze one wavefunction to another.
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
        subroutine orthogonalize_psi (ilogfile, issh, mesh, psi1, psi2, dr)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ilogfile     ! output file unit number
        integer, intent (in) :: issh         ! shell number
        integer, intent (in) :: mesh         ! number of mesh points

        real, intent (in) :: dr              ! grid spacing

        real, intent (in), dimension (mesh) :: psi1
        real, intent (inout), dimension (mesh) :: psi2

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint

        real integral                     ! answer of the psi1*psi2 integration
        real r                            ! value of r at grid point
        real xnorm                        ! normalization of psi2 after

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
        allocate (integrand (mesh))

! Procedure
! ===========================================================================
! First step - integrate the two wavefunctions together.
! The goal is that int(psi1*psi2) = 0 in the end.
        r = - dr
        do ipoint = 1, mesh
          r = r + dr
          integrand(ipoint) = psi1(ipoint)*psi2(ipoint)
        end do
        integral = simpson (mesh, integrand, dr)

! Now orthogonalize the excited state wavefunction (psi2) in a Gram-Schmidt
! orthogonalization algorithm.
        psi2 = psi2 - integral*psi1

! Normalization
        r = - dr
        do ipoint = 1, mesh
          r = r + dr
          integrand(ipoint) = psi2(ipoint)**2
        end do
        xnorm = simpson (mesh, integrand, dr)
        psi2 = psi2/sqrt(xnorm)

! Check normalization
        write (ilogfile,*)
        write (ilogfile,*) ' Inside orthogonalize_psi: Checking normalization [NORM(l) should be 1]'
        r = - dr
        do ipoint = 1, mesh
          r = r + dr
          integrand(ipoint) = psi2(ipoint)**2
        end do
        xnorm = simpson (mesh, integrand, dr)
        write (ilogfile,101) issh, xnorm

! Check orthogonalization
        write (ilogfile,*)
        write (ilogfile,*) ' Inside orthogonalize_psi: Checking orthogonalization [NORM(l) should be 0]'
        r = - dr
        do ipoint = 1, mesh
          r = r + dr
          integrand(ipoint) = psi2(ipoint)*psi1(ipoint)
        end do
        xnorm = simpson (mesh, integrand, dr)
        write (ilogfile,101) issh, xnorm

! Deallocate Arrays
! ===========================================================================
        deallocate (integrand)

! Format Statements
! ===========================================================================
101     format (2x, ' NORM (shell = ', i1, ') = ', f16.12)

        return
        end subroutine orthogonalize_psi
