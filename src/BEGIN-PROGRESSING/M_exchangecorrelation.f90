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

! M_exchangecorrelation.f90
! Program Description
! ===========================================================================
! Program contains the module which handles all calculations involving
! the exchange correlations.
!
! The majority of this code is from:
! Martin Fuchs, FHI der MPG, Berlin, 01-1993
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
module M_exchangecorrelation

! Module Procedures
  contains

! ===========================================================================
! calculate_uxc
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This routine is used as a gateway of sorts to calculate the
! exchangecorrelation potential of the atom. Basically this routine will
! read in the information from M_rcatms and M_atomvar_basic to determine
! which exchangecorrelation is to be calculated, and then pass it along to
! other routines.
!
! The options thus far are:
!
!        1   LDA  Wigner
!        2   LDA  Hedin/Lundqvist
!        3   LDA  Ceperley/Alder Perdew/Zunger (1980)
!        4   GGA  Perdew/Wang (1991)
!        5   GGA  Becke (1988) X, Perdew (1986)
!        6   GGA  Perdew/Burke/Ernzerhof (1996)
!        7   LDA  Zhao/Parr
!        8   LDA  Ceperley/Alder Perdew/Wang (1991)
!        9   GGA  Becke (1988) X, Lee/Yang/Parr (1988)
!        10  GGA  Perdew/Wang (1991) X, Lee/Yang/Parr (1988)
!        11  LSDA Volko/Wilk/Nusair (1980)
!        12  B3LYP  mixing exact exchange and BLYP (9 GGA)
!
! The original code is from
! Martin Fuchs, FHI der MPG, Berlin, 01-1993
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
  subroutine calculate_uxc (mesh, dr, ioption, rho, rhop, rhopp, uxc, exc, ienergy,   &
    &  exmix)
    implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
    integer, intent (in) :: mesh
    real, intent (in) :: dr
    integer, intent (in) :: ienergy
    integer, intent (in) :: ioption
    real, intent (in) :: exmix
    real, intent (in), dimension (mesh) :: rho
    real, intent (in), dimension (mesh) :: rhop
    real, intent (in), dimension (mesh) :: rhopp

! Output
    real, intent(out) :: exc
    real, intent(out), dimension(mesh) :: uxc

! Local Parameters and Data Declaration
! ===========================================================================
! Local Variable Declaration and Description
! ===========================================================================
    integer imesh

    real aln
    real ecp
    real ex
    real fx
    real fxc
    real rh
    real rs
    real rsx
    real zeta

    real, dimension (:), allocatable :: r_mesh
    real, dimension (2) :: cpot
    real, dimension (2) :: d
    real, dimension (:), allocatable :: dec
    real, dimension (:), allocatable :: dex
    real, dimension (:), allocatable :: dexc
    real, dimension (2) :: dp
    real, dimension (2) :: dpp
    real, dimension (2) :: xpot

! Allocate Arrays
! ===========================================================================
    allocate (dec(mesh))
    allocate (dex(mesh))
    allocate (dexc(mesh))
    allocate (r_mesh(mesh))

! Procedure
! ===========================================================================
! Build r_mesh
    do imesh = 1, mesh
      r_mesh(imesh) = dr * (imesh - 1)
    end do

! Select which exchangecorrelation potential is to be calculated
    select case (ioption)

! Wigner
    case (1)
      do imesh = 1, mesh
        rh = rho(imesh)
        call wigner (rh, ex, fx, exc, fxc)
        uxc(imesh) = fxc
        dex(imesh) = ex
        dec(imesh) = exc - ex
      end do

! Hedin-Lundqvist
    case (2)
      dex = 0.d0
      uxc = 0.d0
      do imesh = 1, mesh
        rh = rho(imesh)
        if (rh .ne. 0.0d0) then
          rs = 0.62035049d0*rh**(-0.3333333333333333d0)
          rsx = rs/21.0d0
          aln = dlog(1.0d0 + 1.0d0/rsx)
          ecp = aln + (rsx**3*aln - rsx*rsx) + rsx/2 - 1.0d0/3.0d0
          dex(imesh) = -0.458175d0/rs - 0.0225d0*ecp
          uxc(imesh) = -0.6109d0/rs - 0.0225d0*aln
        end if
      end do

! Ceperley-Alder
! as parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981)
    case (3)
      do imesh = 1, mesh
        rh = rho(imesh)
        call cepal (rh, ex, fx, exc, fxc)
        uxc(imesh) = fxc
        dex(imesh) = ex
        dec(imesh) = exc - ex
      end do

! Perdew exchange, Perdew correlation, generalized gradient approximation 1992
    case (4)
      do imesh = 1, mesh
        d(1) = 0.5d0*rho(imesh)
        dp(1) = 0.5d0*rhop(imesh)
        dpp(1) = 0.5d0*rhopp(imesh)
        d(2) = 0.5d0*rho(imesh)
        dp(2) = 0.5d0*rhop(imesh)
        dpp(2) = 0.5d0*rhopp(imesh)
        call ggaxrad (3, r_mesh(imesh), d, dp, dpp, xpot, dex(imesh))
        call ggacrad (2, r_mesh(imesh), d, dp, dpp, cpot, dec(imesh))
        uxc(imesh) = xpot(1) + cpot(1)
      end do

! Becke exchange, Perdew correlation, generalized gradient approximation
    case (5)
      do imesh = 1, mesh
        d(1) = 0.5d0*rho(imesh)
        dp(1) = 0.5d0*rhop(imesh)
        dpp(1) = 0.5d0*rhopp(imesh)
        d(2) = 0.5d0*rho(imesh)
        dp(2) = 0.5d0*rhop(imesh)
        dpp(2) = 0.5d0*rhopp(imesh)
        call ggaxrad (2, r_mesh(imesh), d, dp, dpp, xpot, dex(imesh))
        call ggacrad (3, r_mesh(imesh), d, dp, dpp, cpot, dec(imesh))
        uxc(imesh) = xpot(1) + cpot(1)
      end do

! Burke-Perdew-Ernzerhof generalized gradient approximation (1996)
    case (6)
      do imesh = 1, mesh
        d(1) = 0.5d0*rho(imesh)
        dp(1) = 0.5d0*rhop(imesh)
        dpp(1) = 0.5d0*rhopp(imesh)
        d(2) = 0.5d0*rho(imesh)
        dp(2) = 0.5d0*rhop(imesh)
        dpp(2) = 0.5d0*rhopp(imesh)
        call ggaxrad (5, r_mesh(imesh), d, dp, dpp, xpot, dex(imesh))
        call ggacrad (5, r_mesh(imesh), d, dp, dpp, cpot, dec(imesh))
        uxc(imesh) = xpot(1) + cpot(1)
      end do

! Wigner-scaled LDA of PRA 46, R5320 (1992)
    case (7)
      do imesh = 1, mesh
        rh = rho(imesh)
        call wigscaled (rh, ex, fx, exc, fxc)
        uxc(imesh) = fxc
        dex(imesh) = ex
        dec(imesh) = exc - ex
      end do

! Ceperley-Alder as parameterized by Perdew-Wang (1991)
    case (8)
      dp(1) = 0.0d0
      dpp(1) = 0.0d0
      dp(2) = 0.0d0
      dpp(2) = 0.0d0
      do imesh = 1, mesh
        d(1) = 0.5d0*rho(imesh)
        d(2) = 0.5d0*rho(imesh)
        call ggaxrad (1, r_mesh(imesh), d, dp, dpp, xpot, dex(imesh))
        call ggacrad (1, r_mesh(imesh), d, dp, dpp, cpot, dec(imesh))
        uxc(imesh) = xpot(1) + cpot(1)
      end do

! Lee-Yang-Parr correlation
! Becke exchange, generalized gradient approximation
    case (9)
      do imesh = 1, mesh
        d(1) = 0.5d0*rho(imesh)
        dp(1) = 0.5d0*rhop(imesh)
        dpp(1) = 0.5d0*rhopp(imesh)
        d(2) = 0.5d0*rho(imesh)
        dp(2) = 0.5d0*rhop(imesh)
        dpp(2) = 0.5d0*rhopp(imesh)
        call ggaxrad (2, r_mesh(imesh), d, dp, dpp, xpot, dex(imesh))
        call ggacrad (4, r_mesh(imesh), d, dp, dpp, cpot, dec(imesh))
        uxc(imesh) = xpot(1) + cpot(1)
      end do

! Perdew-Wang exchange, generalized gradient approximation
    case (10)
      do imesh = 1, mesh
        d(1) = 0.5d0*rho(imesh)
        dp(1) = 0.5d0*rhop(imesh)
        dpp(1) = 0.5d0*rhopp(imesh)
        d(2) = 0.5d0*rho(imesh)
        dp(2) = 0.5d0*rhop(imesh)
        dpp(2) = 0.5d0*rhopp(imesh)
        call ggaxrad (3, r_mesh(imesh), d, dp, dpp, xpot, dex(imesh))
        call ggacrad (4, r_mesh(imesh), d, dp, dpp, cpot, dec(imesh))
        uxc(imesh) = xpot(1) + cpot(1)
      end do

! Volko-Wilk-Nusair exchangecorrelation (spin-polarization)
    case (11)
      do imesh = 1, mesh
        zeta = 0.0d0
        d(1) = 0.5d0*rho(imesh)*(1.0d0 + zeta)
        d(2) = 0.5d0*rho(imesh)*(1.0d0 - zeta)
        call lsdavwn (d, dex(imesh), dec(imesh), xpot, cpot)
        uxc(imesh) = xpot(1) + cpot(1)
      end do

! Lee-Yang-Parr correlation
! Becke exchange, generalized gradient approximation mixed with exact exchange.
    case (12)
      do imesh = 1, mesh
        d(1) = 0.5d0*rho(imesh)
        dp(1) = 0.5d0*rhop(imesh)
        dpp(1) = 0.5d0*rhopp(imesh)
        d(2) = 0.5d0*rho(imesh)
        dp(2) = 0.5d0*rhop(imesh)
        dpp(2) = 0.5d0*rhopp(imesh)
        call ggaxrad (2, r_mesh(imesh), d, dp, dpp, xpot, dex(imesh))
        call ggacrad (4, r_mesh(imesh), d, dp, dpp, cpot, dec(imesh))
        uxc(imesh) = (1.0d0 - exmix)*xpot(1) + cpot(1)
      end do
    end select

! ****************************************************************************
! Take care of the endpoints:
    uxc(1) = 2*uxc(2) - uxc(3)
    uxc(mesh) = 2*uxc(mesh-1) - uxc(mesh - 2)

! Calculate the total energy components
    if (ienergy .eq. 1) then
      if (ioption .ne. 12) then
        dexc = dex + dec
      else
        dexc = (1.0d0 - exmix)*dex + dec
      end if
      exc = 0.0d0
      do imesh = 1, mesh
        exc = exc + dr*r_mesh(imesh)**2*rho(imesh)*(dexc(imesh)-uxc(imesh))
      end do
    end if

! Change the energies and potentials into Rydbergs - multiply by 2!
    exc = 2.0d0*exc
    uxc = 2.0d0*uxc

! Deallocate Arrays
! ===========================================================================
    deallocate (dec)
    deallocate (dex)
    deallocate (dexc)

! Format Statements
! ===========================================================================
    return
  end subroutine calculate_uxc

! End Module
! ===========================================================================
end module M_exchangecorrelation

