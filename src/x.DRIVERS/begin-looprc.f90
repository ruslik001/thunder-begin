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

! begin_looprc.f90
! Program Description
! ===========================================================================
!      This is a driver to loop over many cutoffs for the wavefunctions,
! producing an output file called Z.eig which are the eigenvalues for the
! different cutoffs.
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
        program begin_looprc

! /SYSTEM
        use M_species
        use M_atom_functions
        use M_atomPP_functions

! /BEGIN
        use M_looprc
        use M_rcatms

        implicit none

! Parameters and Data Declaration
! ===========================================================================
       real, parameter :: P_abohr = 0.529177d0   !< Bohr radii to Angstrom

! Variable Declaration and Description
! ===========================================================================
        integer ispecies                     ! looping over species
        integer issh                         ! looping over shells

        real rc                              ! cutoff
        real rcbohr                          ! cutoff in Bohr radii

        real, dimension (:), allocatable :: delta
        real, dimension (:), allocatable :: eigenvalue_rc_max

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
        call read_begin_looprc

! ===========================================================================
! ---------------------------------------------------------------------------
!               G E N E R A T E    W A V E F U N C T I O N S
!                      A N D    P O T E N T I A L S
! ---------------------------------------------------------------------------
! ===========================================================================
! First we evaluate what the eigenvalues are for infinity cutoff.
! Here we set infinity = 15.0d0
! Loop over the number of species
        do ispecies = 1, nspecies
          ! Set the new cutoffs to infinity
          do issh = 1, species(ispecies)%nssh
            species(ispecies)%shell(issh)%rcutoff = rc_infinity(ispecies)
            rcbohr = species(ispecies)%shell(issh)%rcutoff
            species(ispecies)%rcutoff_max = max(rcbohr, species(ispecies)%rcutoff_max)
          end do

! Figure out some stuff in Angstroms
          species(ispecies)%rcutoffA_max = -99.0d0
          do issh = 1, species(ispecies)%nssh
            species(ispecies)%shell(issh)%rcutoffA =                        &
     &        species(ispecies)%shell(issh)%rcutoff*P_abohr
            rc = species(ispecies)%shell(issh)%rcutoffA
            species(ispecies)%rcutoffA_max = max(rc, species(ispecies)%rcutoffA_max)
          end do
        end do

! Initialize the wavefunctions and read pseudopotentials.
! The ispecies loop is inside these subroutines
        call initialize_wf
        call read_vPP

        do ispecies = 1, nspecies
           write(*,*) 'Computing for species', ispecies
! Open the eigenfile for output
          open (unit = 31, file = eigfile(ispecies), status = 'unknown')
          write (31, 101)
          write (31, 100)
          write (31, 101)
          write (31, 102)
          write (31, 101)

          allocate (delta (species(ispecies)%nssh))
          allocate (eigenvalue_rc_max (species(ispecies)%nssh))

          call calculate_rcatm (ispecies)
          eigenvalue_rc_max(:) = wf(ispecies)%shell_data(:)%eigenvalue
          call destroy_rcatm (ispecies)

! Initialize all rcutoffs to rc_min
          do issh = 1, species(ispecies)%nssh
            species(ispecies)%shell(issh)%rcutoff = rc_min(ispecies)
          end do

! Case 1: there is only one shell - s only
! ===========================================================================
          if (species(ispecies)%nssh .eq. 1) then
             write(*,*) 'Case 1: there is only one shell - s only'
! Now start looping over the different cutoffs.
! We start with the initial cutoff provided by the looprc input files.
! We end with the maximum cutoff provided by the looprc input files.
            do while (rc_min(ispecies) .le. rc_max(ispecies))
              species(ispecies)%shell(1)%rcutoff = rc_min(ispecies)
              write (ilogfile,*) ' New cutoffs = ', species(ispecies)%shell(1)%rcutoff

! Figure out some stuff in Angstroms
              species(ispecies)%rcutoff_max = species(ispecies)%shell(1)%rcutoff
              species(ispecies)%shell(1)%rcutoffA = species(ispecies)%shell(1)%rcutoff*P_abohr
              species(ispecies)%rcutoffA_max = species(ispecies)%shell(1)%rcutoffA

! Initialize the wavefunctions
              call initialize_wf

              call calculate_rcatm (ispecies)
              delta(1) = wf(ispecies)%shell_data(1)%eigenvalue - eigenvalue_rc_max(1)
              write (31, 201) species(ispecies)%shell(1)%rcutoff,            &
     &          wf(ispecies)%shell_data(1)%eigenvalue, delta(1)
              call destroy_rcatm (ispecies)

              ! Set the new cutoffs for rc_min
              rc_min(ispecies) = rc_min(ispecies) + 0.05d0
            end do

! Case 2: there are two shells - s and p
! ===========================================================================
          else if (species(ispecies)%nssh .eq. 2) then
             write(*,*) 'Case 2: there are two shells - s and p'
! Now start looping over the different cutoffs.
! We start with the initial cutoff provided by the looprc input files.
! We end with the maximum cutoff provided by the looprc input files.
            do while (rc_min(ispecies) .le. rc_max(ispecies))
               write(*,*) 'rc_min', rc_min(ispecies)
              ! Set up initial cutoffs
              do issh = 1, species(ispecies)%nssh
                species(ispecies)%shell(issh)%rcutoff = rc_min(ispecies)
              end do

              ! looping over the p states
              do while (species(ispecies)%shell(2)%rcutoff .le. rc_min(ispecies) + 0.8d0)
                 write(*,*) 'Species=', ispecies,' Rcutoff=', species(ispecies)%shell(2)%rcutoff
                write (ilogfile,*) ' New cutoffs = ',                        &
     &            species(ispecies)%shell(1)%rcutoff,                        &
     &            species(ispecies)%shell(2)%rcutoff

! Find the maximum cutoff
                do issh = 1, species(ispecies)%nssh
                   write(*,*) issh
                  rcbohr = species(ispecies)%shell(issh)%rcutoff
                  species(ispecies)%rcutoff_max = max(rcbohr, species(ispecies)%rcutoff_max)
                end do

! Figure out some stuff in Angstroms
                species(ispecies)%rcutoffA_max = -99.0d0
                do issh = 1, species(ispecies)%nssh
                  species(ispecies)%shell(issh)%rcutoffA =                   &
     &            species(ispecies)%shell(issh)%rcutoff*P_abohr
                  rc = species(ispecies)%shell(issh)%rcutoffA
                  species(ispecies)%rcutoffA_max = max(rc, species(ispecies)%rcutoffA_max)
                end do

! Initialize the wavefunctions
                call initialize_wf

                call calculate_rcatm (ispecies)
                do issh = 1, species(ispecies)%nssh
                  delta(issh) = wf(ispecies)%shell_data(issh)%eigenvalue     &
     &                         - eigenvalue_rc_max(issh)
                end do

                write (31, 202) species(ispecies)%shell(1)%rcutoff,          &
     &            wf(ispecies)%shell_data(1)%eigenvalue, delta(1),           &
     &            species(ispecies)%shell(2)%rcutoff,                        &
     &            wf(ispecies)%shell_data(2)%eigenvalue, delta(2)
                call destroy_rcatm (ispecies)

                ! Set new cutoffs for the second shell (p state)
                species(ispecies)%shell(2)%rcutoff =                         &
     &            species(ispecies)%shell(2)%rcutoff + 0.050d0
              end do
              ! Set the new cutoffs for rc_min
              rc_min(ispecies) = rc_min(ispecies) + 0.050d0
            end do

! Case 3: there are three shells - s, p, and d
! ===========================================================================
          else if (species(ispecies)%nssh .eq. 3) then

! Now start looping over the different cutoffs.
! We start with the initial cutoff provided by the looprc input files.
! We end with the maximum cutoff provided by the looprc input files.
            do while (rc_min(ispecies) .le. rc_max(ispecies))

              ! Set up initial cutoffs
              do issh = 1, species(ispecies)%nssh
                species(ispecies)%shell(issh)%rcutoff = rc_min(ispecies)
              end do

              ! looping over the p states
              do while (species(ispecies)%shell(2)%rcutoff .le. rc_min(ispecies) + 0.8d0)
                 write(*,*) 'Loop', species(ispecies)%shell(2)%rcutoff

                ! looping over the d states
                species(ispecies)%shell(3)%rcutoff = species(ispecies)%shell(3)%rcutoff - 0.8d0
                do while (species(ispecies)%shell(3)%rcutoff .le. rc_min(ispecies) + 0.8d0)
                  write (ilogfile,*) ' New cutoffs = ',                      &
     &              species(ispecies)%shell(1)%rcutoff,                      &
     &              species(ispecies)%shell(2)%rcutoff,                      &
     &              species(ispecies)%shell(3)%rcutoff

! Find the maximum cutoff
                  do issh = 1, species(ispecies)%nssh
                    rcbohr = species(ispecies)%shell(issh)%rcutoff
                    species(ispecies)%rcutoff_max = max(rcbohr, species(ispecies)%rcutoff_max)
                  end do

! Figure out some stuff in Angstroms
                  species(ispecies)%rcutoffA_max = -99.0d0
                  do issh = 1, species(ispecies)%nssh
                    species(ispecies)%shell(issh)%rcutoffA =                 &
     &              species(ispecies)%shell(issh)%rcutoff*P_abohr
                    rc = species(ispecies)%shell(issh)%rcutoffA
                    species(ispecies)%rcutoffA_max = max(rc, species(ispecies)%rcutoffA_max)
                  end do

! Initialize the wavefunctions
                  call initialize_wf

                  call calculate_rcatm (ispecies)
                  do issh = 1, species(ispecies)%nssh
                    delta(issh) = wf(ispecies)%shell_data(issh)%eigenvalue   &
     &                           - eigenvalue_rc_max(issh)
                  end do

                  write (31, 203) species(ispecies)%shell(1)%rcutoff,        &
     &              wf(ispecies)%shell_data(1)%eigenvalue, delta(1),         &
     &              species(ispecies)%shell(2)%rcutoff,                      &
     &              wf(ispecies)%shell_data(2)%eigenvalue, delta(2),         &
     &              species(ispecies)%shell(3)%rcutoff,                      &
     &              wf(ispecies)%shell_data(3)%eigenvalue, delta(3)
                  call destroy_rcatm (ispecies)

                  ! Set new cutoffs for the third shell (d state)
                  species(ispecies)%shell(3)%rcutoff =                       &
     &              species(ispecies)%shell(3)%rcutoff + 0.050d0
                end do
                ! Set new cutoffs for the second shell (p state)
                species(ispecies)%shell(2)%rcutoff =                         &
     &            species(ispecies)%shell(2)%rcutoff + 0.050d0
              end do
              ! Set the new cutoffs for rc_min
              rc_min(ispecies) = rc_min(ispecies) + 0.050d0
            end do
          end if ! end different cases

          deallocate (delta, eigenvalue_rc_max)
        end do  ! end loop over species

        call cpu_time (time_end)
        write (ilogfile,*)
        write (ilogfile,*) ' BEGIN RUNTIME : ', time_end - time_begin, '[sec] '
        close (ilogfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (10x, 'SHELL 1', 10x, '|', 10x, 'SHELL 2', 10x, '|', 10x, 'SHELL 3')
101     format (27('='), '+', 27('='), '+', 27('='))
102     format (2('   rc   eigenval    delta  |'), '   rc   eigenval    delta')

201     format (1(f6.3, 2(1x, f9.5), ' |'))
202     format (2(f6.3, 2(1x, f9.5), ' |'))
203     format (2(f6.3, 2(1x, f9.5), ' |'), f6.3, 2(1x, f9.5))

! End Program
! ===========================================================================
        stop
        end program begin_looprc
