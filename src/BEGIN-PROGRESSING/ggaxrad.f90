! copyright info:
!
!                             @Copyright 2008
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

! ggaxrad.f90
! Program Description
! ===========================================================================
!
!      This routine calculates the exchange potential and energy density.
! Spherical symmetry is used. LSDA - GGA
!
! input
!    mode = 1    LSDA
!    mode = 2    GGA-X Becke
!    mode = 3    GGA-X Perdew
!    mode = 5    GGA-X Burke-Perdew-Ernzerhof
!
! convention
!    yy(1) = spin up, yy(2) = spin down
!
! Martin Fuchs, FHI der MPG, Berlin, 07-1992
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
        subroutine ggaxrad (mode, r_in, rho_in, rhop_in, rhopp_in, xpot, xen)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mode
        real, intent (in) :: r_in
        real, intent (in), dimension (2) :: rho_in
        real, intent (in), dimension (2) :: rhop_in
        real, intent (in), dimension (2) :: rhopp_in
! Output
        real, intent (out) :: xen
        real, intent (out), dimension (2) :: xpot

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: eps = 1.0d-15
        real, parameter :: pi = 3.141592653589793238462643d0

! Local Variable Declaration and Description
! ===========================================================================
        integer ispin
        real density
        real densityp
        real densitypp
        real ex
        real fermik
        real r
        real s
        real u
        real v
        real vx

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! If r is really small, then set to manageably small number.
        r = r_in
        if (r_in .lt. 1.0d-4) r = 1.0d-4

! exchange GGA, loop for up & down spin
        xen = 0.0d0
        do ispin = 1, 2
         if (rho_in(ispin) .le. eps) then
          xpot(ispin) = 0.0d0
         else
          density = 2.0d0*rho_in(ispin)
          if (mode .eq. 1) then
           call xlda (density, vx, ex)
          else if (mode .eq. 2 .or. mode .eq. 3 .or. mode .eq. 5) then
           densityp = 2.0d0*rhop_in(ispin)
           densitypp = 2.0d0*rhopp_in(ispin)
           fermik = 2.0d0*(3.0d0*pi*pi*density)**(1.0d0/3.0d0)
           s = abs(densityp)/(fermik*density)
           u = abs(densityp)*densitypp/(density*density*fermik**3)
           v = (densitypp + 2.0d0*densityp/r)/(density*fermik*fermik)
           select case (mode)
            case (2)
             call xbecke (density, s, u, v, ex, vx)
            case (3)
             call exch (density, s, u, v, ex, vx)
            case (5)
             call exchpbe (density, s, u, v, 1, 1, ex, vx)
           end select
          else
           stop 'ggaxrad : mode improper'
          end if
          xpot(ispin) = vx
          xen = xen + rho_in(ispin)*ex
         end if
        end do

! Energy
        xen = xen/max(rho_in(1) + rho_in(2), eps)

! Deallocate Arrays
! ===========================================================================
! Format Statements
! ===========================================================================

        return
        end subroutine ggaxrad
