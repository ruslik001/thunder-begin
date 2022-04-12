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

! ggacrad.f90
! Program Description
! ===========================================================================
!
!      This routine calculates the correlation potential and energy density.
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
        subroutine ggacrad (mode, r_in, rho_in, rhop_in, rhopp_in, cpot, cen)
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
        real, intent (out) :: cen
        real, intent (out), dimension (2) :: cpot

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: pi = 3.141592653589793238462643d0
        real, parameter :: crs = 1.91915829267751281d0
        real, parameter :: eps = 1.0d-15

! Local Variable Declaration and Description
! ===========================================================================
        real alfc
        real density
        real densityp
        real densityp11
        real densityp12
        real densityp22
        real densitypp
        real ec
        real ecrs
        real eczet
        real fermik
        real g
        real gsfermik
        real h
        real r
        real rs
        real sfermik
        real t
        real uu
        real vv
        real ww
        real zet
        real ztp

        real, dimension (2) :: dvc, vc
        real, dimension (2) :: flip

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! If r is really small, then set to manageably small number.
        r = r_in
        if (r_in .lt. 1.0d-4) r = 1.0d-4

! LSDA
        density = rho_in(1) + rho_in(2)
        densityp = rhop_in(1) + rhop_in(2)
        densitypp = rhopp_in(1) + rhopp_in(2)
        cen = 0.0d0
        if (density .le. eps) then
         cen = 0.0d0
         cpot(1) = 0.0d0
         cpot(2) = 0.0d0
         return
        end if

        if (mode .ne. 4) then
         zet = (rho_in(1) - rho_in(2))/density
         ztp = (rhop_in(1) - rhop_in(2) - zet*densityp)/density
         fermik = (3.0d0*pi*pi*density)**(1.0d0/3.0d0)
         rs = crs/fermik
         call corlsd (rs, zet, ec, vc(1), vc(2), ecrs, eczet, alfc)

! GGA correction to LSDA
         select case (mode)
          case (1)
           dvc = 0.0d0
           h = 0.0d0
          case (2)
           sfermik = 2.0d0*sqrt(fermik/pi)
           g = ((1.0d0 + zet)**(2.0d0/3.0d0)                               &
     &        + (1.0d0 - zet)**(2.0d0/3.0d0))/2.0d0
           gsfermik = 2.0d0*sfermik*g
           t = abs(densityp)/(density*gsfermik)
           uu = abs(densityp)*densitypp/(density*density*gsfermik**3)
           vv = (densitypp + 2.0d0*densityp/r)/(density*gsfermik*gsfermik)
           ww = densityp*ztp/(density*gsfermik*gsfermik)
           call corgga (rs, zet, t, uu, vv, ww, h, dvc(1), dvc(2))
          case (3)
           uu = abs(densityp)*densitypp
           vv = densitypp + 2.0d0*densityp/r
           densityp11 = rhop_in(1)*rhop_in(1)
           densityp22 = rhop_in(2)*rhop_in(2)
           densityp12 = rhop_in(1)*rhop_in(2)
           if (rhop_in(1) .ne. 0.0d0 .or. rhop_in(2) .ne. 0.0d0)then
            call corga86 (rho_in(1),rho_in(2), densityp11, densityp22,     &
     &                    densityp12, uu, vv, h, dvc(1), dvc(2))
           end if
          case (5)
           sfermik = 2.0d0*sqrt(fermik/pi)
           g = ((1.0d0 + zet)**(2.0d0/3.0d0)                               &
     &        + (1.0d0 - zet)**(2.0d0/3.0d0))/2.0d0
           gsfermik = 2.0d0*sfermik*g
           t = abs(densityp)/(density*gsfermik)
           uu = abs(densityp)*densitypp/(density*density*gsfermik**3)
           vv = (densitypp + 2.0d0*densityp/r)/(density*gsfermik*gsfermik)
           ww = densityp*ztp/(density*gsfermik*gsfermik)
           call corpbe (rs, zet, t, uu, vv, ww, 1, 1, ec, vc(1), vc(2), h,  &
     &                  dvc(1),dvc(2))
          end select
          cpot = vc + dvc
          cen = ec + h
        else if (mode .eq. 4) then
         flip = rhopp_in + 2.0d0*rhop_in/r
         call corlyp (.true., rho_in(1), rho_in(2), rhop_in(1), rhop_in(2), &
     &                flip(1), flip(2), cen, cpot(1), cpot(2))
        else
         stop 'ggacrad : mode improper'
        end if

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end subroutine ggacrad

