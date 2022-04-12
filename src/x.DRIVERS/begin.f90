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

! begin.f90
! Program Description
! ===========================================================================
!      Driver to run FIREBALL version which includes BEGIN and CREATE.
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
        program begin

! /SYSTEM
        use M_species
        use M_atom_functions
        use M_atomPP_functions

! /BEGIN
        use M_rcatms
        use M_rcatms_excited
        use M_vnn

        implicit none

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies                     ! looping over species

! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_begin
        real time_end

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! ===========================================================================
! ---------------------------------------------------------------------------
!                             W E L C O M E !
! ---------------------------------------------------------------------------
! ===========================================================================
        call cpu_time (time_begin)
        open (unit = ilogfile, file = 'output.log', status = 'replace')
        call welcome_begin

! ===========================================================================
! ---------------------------------------------------------------------------
!             R E A D   I N   S Y S T E M   I N F O R M A T I O N
! ---------------------------------------------------------------------------
! ===========================================================================
! Call read_input to define variables from the begin.input file
        call read_Fdata_location
        allocate (species_PP (nspecies))
        call read_begin
        call write_create

! ===========================================================================
! ---------------------------------------------------------------------------
!               G E N E R A T E    W A V E F U N C T I O N S
! ---------------------------------------------------------------------------
! ===========================================================================
! Initialize the wavefunctions and read pseudopotentials.
        call initialize_wf
        call read_vPP

! Read additional PP files and initialize excited state wavefunctions.
        call initialize_wf_ion
        call read_vPP_ion

! Loop over the number of species
        do ispecies = 1, nspecies

! nexcite = 0: find the ground state only
          if (species(ispecies)%nexcite .eq. 0) then
            call calculate_rcatm (ispecies)
            call writeout_wf (ispecies)
            call destroy_rcatm (ispecies)

! nexcite = 1: find the ground state, psi with DMOL formalism, orthogonalize
          else if (species(ispecies)%nexcite .eq. 1) then
            call calculate_rcatm (ispecies)
            call writeout_wf (ispecies)
            call calculate_rcatm_excited (ispecies)

            call destroy_rcatm (ispecies)
            call destroy_rcatm_excited (ispecies)
          end if
        end do

! ===========================================================================
! ---------------------------------------------------------------------------
!               G E N E R A T E    P O T E N T I A L S
! ---------------------------------------------------------------------------
! ===========================================================================
        call read_create
        call read_wavefunctions
        call initialize_na

! Loop over the number of species
        do ispecies = 1, nspecies
          call calculate_vnn (ispecies)
          call destroy_na (ispecies)
        end do

        call cpu_time (time_end)
        write (ilogfile,*)
        write (ilogfile,*) ' BEGIN RUNTIME : ', time_end - time_begin, '[sec] '
        close (ilogfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Program
! ===========================================================================
        stop
        end program begin
