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

! corlyp.f90
! Program Description
! ===========================================================================
!       Lee Yang Parr correlation energy functional one-dimensional densities
! only no provisions taken against division by zero.
!
! See e.g. C. Lee et al. Phys. Rev. B 37 785 (1988)
!
! Hartree A.U.
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
! Program Declaration
! ===========================================================================
        subroutine corlyp (tpot, dp, dm, dp1, dm1, dp2, dm2, ec, vcp0, vcm0)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: dm, dm1, dm2
        real, intent (in) :: dp, dp1, dp2

        logical, intent (in) :: tpot

! Output
        real, intent (out) :: ec
        real, intent (out) :: vcm0
        real, intent (out) :: vcp0

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: aa = 0.04918d0
        real, parameter :: bb = 0.132d0
        real, parameter :: cc = 0.2533d0
        real, parameter :: dd = 0.349d0
        real, parameter :: c5 = 4.55779986d0
        real, parameter :: c6 = 1.0d0/72.0d0
        real, parameter :: c7 = 1.0d0/18.0d0
        real, parameter :: c8 = 0.125d0
        real, parameter :: t13 = 1.0d0/3.0d0
        real, parameter :: t89 = 8.0d0/9.0d0

! Local Variable Declaration and Description
! ===========================================================================
        real c1, c2, c3, c4, c9
        real chf
        real d0xt13, d0xt53
        real d0, d1, d2
        real dmt53, dpt53
        real dxsq
        real ga
        real gafm, gafp
        real gb
        real h
        real h2
        real hf
        real hff
        real sc
        real sc2
        real scf
        real t43, t53, t83
        real yafm, yafp
        real yb, yb1, yb2
        real ybfm, ybfp
        real yy1
        real yz, yz1, yz2
        real z1, z2
        real zfm, zfp
        real zeta

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Initialize some parameters
        c1 = -4.0d0*aa
        c2 = dd
        c3 = 2.0d0*bb
        c4 = cc
        t53 = 5.0d0*t13
        t43 = 4.0d0*t13
        t83 = 2.0d0*t43
        c9 = t43 + t89

        d0 = dp + dm
        dxsq = 1.0d0/(d0*d0)
        d1 = dp1 + dm1
        d2 = dp2 + dm2
        d0xt13 = d0**(-t13)
        d0xt53 = d0xt13*d0xt13/d0
        dpt53 = dp**t53
        dmt53 = dm**t53

! Polarization factor
        zeta = c1*(dp*dm)*dxsq

! Scaling function
        sc = 1.0d0/(1.0d0 + c2*d0xt13)
        h = c3*d0xt53*exp(-c4*d0xt13)

! kinetic energy density expansion
        ga = c5*(dp*dpt53 + dm*dmt53)
        gb = c6*(dp1*dp1 - dp*dp2 + dm1*dm1 - dm*dm2) + c7*(dp*dp2 + dm*dm2) &
     &      + c8*(d0*d2 - d1*d1)

! Calculate potential
        if (tpot) then
         gafp = t83*c5*dpt53
         gafm = t83*c5*dmt53

         scf = t13*c2*d0xt13/d0*sc*sc
         sc2 = scf*(d2 + 2.0d0*(scf/sc - 2.0d0*t13/d0)*d1*d1)

         chf = t13*(c4*d0xt13 - 5.0d0)/d0
         hf = chf*h
         hff = h*(chf**2 + t13*(5.0d0 - 4.0d0*t13*c4*d0xt13)*dxsq)
         h2 = (hf*d2 + hff*d1*d1)

         zfp = (c1*dm - 2.0d0*zeta*d0)*dxsq
         zfm = (c1*dp - 2.0d0*zeta*d0)*dxsq
         yz = zeta/c1
         yy1 = dp*dm1 + dm*dp1
         yz1 = (yy1 - 2.0d0*yz*d1*d0)*dxsq
         yz2 = (2.0d0*yz*d1*d1 - 2.0d0*(yz1*d1 + yz*d2)*d0                   &
     &        - 2.0d0*d1*yy1/d0 + (dp*dm2 + 2.0d0*dp1*dm1 + dm*dp2))*dxsq
         z1 = c1*yz1
         z2 = c1*yz2

         yafp = sc*(d0*zfp + zeta) + zeta*d0*scf
         yafm = sc*(d0*zfm + zeta) + zeta*d0*scf

         yb = sc*zeta*h
         ybfp = sc*(h*zfp + zeta*hf) + zeta*h*scf
         ybfm = sc*(h*zfm + zeta*hf) + zeta*h*scf
         yb1 = sc*(h*z1 + zeta*hf*d1) + zeta*h*scf*d1
         yb2 = (sc*hf + h*scf)*d1*z1 + h*sc*z2 + (sc*z1 + zeta*scf*d1)*hf*d1 &
     &        + zeta*sc*h2 + (zeta*hf*d1 + h*z1)*scf*d1 + zeta*h*sc2

! Collect contributions
         vcp0 = yafp + ybfp*(ga + gb)                                        &
     &         + yb*(gafp + 2.0d0*c8*(c9*dp2 + 2.0d0*dm2))                   &
     &         + yb1*2.0d0*c8*(c9*dp1 + 2.0d0*dm1) + yb2*c8*(t43*dp + dm)
         vcm0 = yafm + ybfm*(ga + gb)                                        &
     &         + yb*(gafm + 2.0d0*c8*(c9*dm2 + 2.0d0*dp2))                   &
     &         + yb1*2.0d0*c8*(c9*dm1 + 2.0d0*dp1) + yb2*c8*(t43*dm + dp)
        else
         vcp0 = 0.0d0
         vcm0 = 0.0d0
        endif

! Correlation energy per electron
        ec = zeta*sc*(d0 + h*(ga + gb))/d0

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end
