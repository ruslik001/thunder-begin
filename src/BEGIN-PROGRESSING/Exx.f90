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
!
!     ***********************************************
!
      SUBROUTINE ExxPot(mx,ms,mmax,dr,rp,Rall,l,ir,i1,i2,occ,cgcoeff,vexx)
!
!
!     ***********************************************
!
! Argument Declaration and Description
! ===========================================================================
! Input
!
      integer, intent (in) ::  mx
      integer, intent (in) ::  ms
      integer, intent (in) ::  mmax
      integer, intent (in) ::  ir
      integer, intent (in) ::  i1
      integer, intent (in) ::  i2

      integer, intent (in), dimension (ms) :: l

      real :: dr
      real, intent (in), dimension (0:3,0:3,0:6) :: cgcoeff
      real, intent (in), dimension (mx) :: rp
      real, intent (in), dimension (ms) :: occ
      real, intent (in), dimension (ms,mx) :: Rall

! Output
!
      real, intent (out), dimension (ms,mx) :: vexx
!
!
! Local Parameters and Data Declaration
! =========================================================================
!
!
      integer :: max,i,il,k
      real, dimension(:), allocatable :: rh
      real  :: f1small,f2small
      real  :: f1large,f2large

! =========================================================================
!
!     This routine computes the exact exchange potential
!     for the state Rall(n(ir),l(ir)).
!     All quantities are given on the linear mesh
!     points rp(i).
!
!     INPUT:
!
!       mmax        size of the fields
!       dr          r(2) - r(1)
!       rp(*)       radial mesh
!       Rall(*,*)   all wave functions  ( = uall(*,i)/rp(i) )
!       l(*)        angular quantum number
!       ir          index of orbital, for which vxch is computed
!       i1          index of first occupied state in sc loop
!       i2          index of highest occupied state
!       occ(*)      number of electrons in shell
!
!       cgcoeff(ll,lr,k)  Clebsch Gordan coefficients
!                         < ll,0  k,0  |  lr,0 >**2
!
!     OUTPUT:
!
!       vexx(*,*) exchange potential for
!                 all states
!
!     **********************************************************
!
!       The integral has the form
!
!     int dr' [ I(r,r') * r'*r']  (from r(1) to r(mmax) )
!
!     I(r,r') =
!     -1.0*sum_j occ(j) Rall(j,r') * Rall(ir,r') * Rall(j,r) *
!          sum_k < ll,0  k,0  |  lr,0 >**2 ( r_< / r_> )^k / r_>
!
!     The potential is the integral devided by
!     Rall(ir,r)*2(2*l(ir)+1)
!
!
!     The cases r_<  and r_> lead to a two regions
!     integral. This can be solved by two loops, one
!     over i=1,mmax, the other over i=mmax,1,-1
!     Integration by trapezoid rule
!
!
!     *****************************************
!     Written by Juergen Fritsch
!
      allocate (rh(mx))  !  the field has to be deallocated at the end
!
      rh(:) = rp(:)
      rh(1) = 1.0E-7
!
!
      max=mmax
      if(mod(mmax,2) .eq. 0) then
        write(*,*) ' vExch - got even no. of mesh points max, assume mmax-1'
        max=max-1
      endif
!
      DO i=1,max
        vexx(ir,i) = 0.0
      END DO
!
      DO il=i1,i2    ! :::::::::::::  il loop  ::::::::::
        DO k=0,l(il)+l(ir)  !  -----   k loop   ---------
!
          f1small = 0.0
          DO i=2,max
!
            f2small = 0.5*dr* &
              & occ(il)*Rall(il,i)*Rall(ir,i)* &
              & cgcoeff(l(il),l(ir),k)*rh(i)**(k+2)
!
            f1small = f1small + f2small
            vexx(ir,i) = vexx(ir,i) + &
              & f1small*Rall(il,i)/(rh(i)**(k+1))
            f1small = f1small + f2small
          END DO
!
          f1large = 0.0
          DO i=max,2,-1
!
            f2large = 0.5*dr* &
              & occ(il)*Rall(il,i)*Rall(ir,i)* &
              & cgcoeff(l(il),l(ir),k)*(rh(i)**2/rh(i)**(k+1))
!
            f1large = f1large + f2large
            vexx(ir,i) = vexx(ir,i) +  &
              & f1large*Rall(il,i)*(rh(i)**k)
            f1large = f1large + f2large
!
          END DO
!
!         -----------------
!         Finish the integration up to r(1)
!
          vexx(ir,1) = vexx(ir,1) + &
                     & f1large*Rall(il,1)*(rh(1)**k)
!         -----------------
!
!
        END DO  !  --------  k  loop  ------------------
      END DO    !  ::::::::: il lool  ::::::::::::::::::
!
!
!     Finalize
!
      DO i=1,max
         IF ( abs(Rall(ir,i)) .GT. 1.0E-4) THEN
           vexx(ir,i) = -2.0*vexx(ir,i) / ( Rall(ir,i)*2.0*(2.0*l(ir)+1.0) )
         ELSE
           vexx(ir,i) = 0.0
         END IF
!
      END DO
!
!     Factor 2.0 is for conversion to Rydberg (from Hartree)
!
      deallocate(rh)
!
      RETURN
      END
