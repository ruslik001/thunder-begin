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

! M_looprc
! Module Description
! ===========================================================================
!>       This is a module that contains information for the looprc program.
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
        module M_looprc
        use M_species

! Type Declaration
! ===========================================================================
! For looprc information
        real, dimension (:), allocatable :: rc_min   ! looprc minimum cutoff
        real, dimension (:), allocatable :: rc_max   ! looprc maximum cutoff
        real, dimension (:), allocatable :: rc_infinity  ! infinity cutoff

        ! rcutoff eigenvalues file
        character (len = 7), dimension (:), allocatable :: eigfile

! Variable Declaration and Description
! ===========================================================================
! None

! module procedures
        contains


! ===========================================================================
! read_begin_looprc
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine will read the begin.input file and define the
! variables assoiated with the type of species and exchange correlation.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine read_begin_looprc
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
        character*2 periodic (103)
        data periodic  / 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ',    &
             'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar',    &
             'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',    &
             'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',    &
             'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',    &
             'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce',    &
             'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',    &
             'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt',    &
             'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra',    &
             'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',    &
             'Es', 'Fm', 'Md', 'No', 'Lw' /

! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint                     ! counter over character point
        integer iten, ione                 ! character array places
        integer ispecies                   ! counter over species
        integer issh                       ! counter over shells
        integer lssh                              ! l quantum number of shell
        integer nssh                       ! number of shells
        integer nzx_max                    ! what is the maximum Z
        integer nZ                         ! atomic number

        logical read_input
        logical file_exists

        character (len = 2) atomcheck
        character (len = 11) buffer        ! buffer for generating wavefunction
        character (len = 15) inputfile     ! input file for species info

        character (len = 1), dimension (0:9) :: z

! Allocate Arrays
! ===========================================================================
        allocate (rc_min (nspecies))
        allocate (rc_max (nspecies))
        allocate (rc_infinity (nspecies))
        allocate (eigfile (nspecies))

! Procedure
! ===========================================================================
! Initialize some symbols
        do ipoint = 0, 9
          z(ipoint) = char(48 + ipoint)
        end do

! We now read in a begin.input file. This determines the number of atoms
! and the types of atoms.
        write (ilogfile,*) ' We now read looprc.inp '
        INQUIRE(FILE="looprc.inp", EXIST=file_exists)   ! file_exists will be TRUE if the file
                                                        ! exists and FALSE otherwise
        if ( file_exists ) then
           open (unit = 11, file = 'looprc.inp', status = 'old')
        else
           write(*,*) 'ERROR: Could not open: "looprc.inp"'
           call exit(1)
        end if

        read (11, 101) signature

        nzx_max = 0
        do ispecies = 1, nspecies
          read (11,101) inputfile
          inquire (file = inputfile, exist = read_input)
          if (read_input) then
            open (unit = 12, file = inputfile, status = 'old')
          else
            write (ilogfile,*) ' The following input file does not exist! '
            write (ilogfile,102) inputfile
            stop ' Error in read_begin_looprc! '
          end if

          read (12,101) species(ispecies)%name
          read (12,103) species(ispecies)%symbol
          read (12,*) nZ
          if (nZ .ne. species(ispecies)%nZ) then
            stop ' inconsistency between Fdata.inp and create.inp in nZ '
          end if
          nzx_max = max(species(ispecies)%nZ,nzx_max)
          if (species(ispecies)%nZ .lt. nzx_max) then
            write (ilogfile,*) ' ispecies = ', ispecies
            write (ilogfile,*) ' Z(ispecies) .lt. Z(ispecies-1) '
            write (ilogfile,*) ' nzx(ispecies) = ', species(ispecies)%nZ
            write (ilogfile,*) ' nzx(ispecies-1) = ', species(ispecies-1)%nZ
            stop ' Must stop bad order.  Z1 < Z2 < Z3 ... violated'
          end if

! Check whether you put in the correct nz for that atom.
! Go through the periodic table and check.
          atomcheck = periodic(species(ispecies)%nZ)
          write (ilogfile,104) ispecies, species(ispecies)%nZ,               &
     &                 species(ispecies)%symbol, atomcheck
          if (species(ispecies)%symbol .ne. atomcheck) stop                  &
     &      ' Wrong nZ(nuc) for atom!!'

! The variable symbol is a three character string which has the atomic number.
! For example, 001 is hydrogen, 008 is oxygen, and 094 is plutonium.
          nZ = species(ispecies)%nZ
          iten = nZ/10
          if (nZ .lt. 10) buffer(1:3) = z(0)//z(0)//z(nZ)
          if (nZ .ge. 10 .and. nZ .le. 99) then
            iten = nZ/10
            ione = nZ - 10*iten
            buffer(1:3) = z(0)//z(iten)//z(ione)
          end if

! Read mass
          read (12,*) species(ispecies)%xmass

! Read number of valence electrons
          read (12,*) species(ispecies)%Zval

! For the pseudopotential filename
          species(ispecies)%PPfile(1:6) = 'basis/'
          species(ispecies)%PPfile(7:9) = buffer(1:3)
          species(ispecies)%PPfile(10:12) = '.pp'

! For the eigfile output
          eigfile(ispecies)(1:3) = buffer(1:3)
          eigfile(ispecies)(4:7) = '.eig'

! Read in the exchange-correlation (which is a number)
          read (12,*) species(ispecies)%ixc_option

! Now read stuff in related to the orbital information
          read (12,*) nssh
          species(ispecies)%nssh = nssh
          allocate (species(ispecies)%shell(nssh))

! Loop over the number of shells:
! Read the l quantum number (0, 1, 2, and 3 => s, p, d, and f), the occupation
! number, cutoff radius (in bohr), wavefunction, and neutral atom potential
! for each shell.
          species(ispecies)%rcutoff_max = -99.0d0
          do issh = 1, species(ispecies)%nssh
            read (12,*) species(ispecies)%shell(issh)%lssh
            read (12,*) species(ispecies)%shell(issh)%Qneutral

! What is a0?  It is the initial guess for the wavefunctions.
! We assume a hydrogenic atom, so psi = exp(-a0/r)
! Set the values here.
            lssh = species(ispecies)%shell(issh)%lssh
            if (lssh .eq. 0) species(ispecies)%shell(issh)%a0 = 2.0d0
            if (lssh .eq. 1) species(ispecies)%shell(issh)%a0 = 1.0d0
            if (lssh .eq. 2) species(ispecies)%shell(issh)%a0 = 0.8d0
          end do

          read (12,*) rc_min(ispecies)
          read (12,*) rc_max(ispecies)
          read (12,*) rc_infinity(ispecies)

! End loop over species
          close (unit = 12)
        end do
        close (unit = 11)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (a15)
102     format (2x, a15)
103     format (a2)
104     format (2x, ' Species = ', i2, ' Nuclear Z = ', i3, ' atomname = ',  &
     &           a2, ' in periodic table = ', a2)

! ===========================================================================
        return
        end subroutine read_begin_looprc


! End Module
! ===========================================================================
        end module M_looprc
