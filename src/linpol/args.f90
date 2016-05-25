!#########################################
!# LINPOL: LINear POLymer instanton code #
!#########################################
!#
!# Research conducted using LINPOL or any derivative thereof should cite the following PhD thesis:
!#
!# A. A. Reid, "Quantum Tunnelling Splittings in Water Clusters, from Ring-Polymer Instanton Theory", University of Cambridge (2014)
!#
!#########################################
!#
!# This file is part of LINPOL, Copyright (C) Adam Reid, Jeremy Richardson & Michael F. Herbst 2011-2016.
!#
!# LINPOL is free software: you can redistribute it and/or modify
!# it under the terms of the GNU General Public License as published by
!# the Free Software Foundation, either version 3 of the License, or
!# (at your option) any later version.
!#
!# LINPOL is distributed in the hope that it will be useful,
!# but WITHOUT ANY WARRANTY; without even the implied warranty of
!# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!# GNU General Public License for more details.
!#
!# You should have received a copy of the GNU General Public License
!# along with LINPOL.  If not, see <http://www.gnu.org/licenses/>.
!#
!#########################################
module args
use type_module
use parameters
use pes_interface
use opt_interface
use opt_selected
implicit none
private         !everything private
!
! =======================================================================
! PURPOSE OF THIS MODULE
! =======================================================================
! Parse command line args and save to parameters 
!
! Record of revisions:
!
!  Date           Person         Change implemented
!  ----           ------         ------------------
!  20 Nov 2012    M. Herbst      Development of module
!  16 Jul 2013    A. Reid        Addition of hesseps option
!
! =======================================================================

public :: parse_args                    ! The only public subroutine

integer :: myrank, error                ! For MPI
integer :: argc                         ! No. of command line arguments

contains
subroutine print_help(verbous)
        logical,intent(in) :: verbous
        CHARACTER(len=10) :: loosedefault,fixeddefault
        if (myrank .ne. root) return

        if (endchoice .eq. 2) then
                loosedefault=" (DEFAULT)"
                fixeddefault=""
        else
                loosedefault=""
                fixeddefault=" (DEFAULT)"
        end if

        write(*,*) 'linpolexe [ Options ] <INPUT>'
        write(*,*) 'Program to optimise an initial guess for a linear polymer instanton.'
        write(*,*) 'INPUT is a valid xyz file with the starting geometry.'
        write(*,*) ''
        write(*,*) 'Options:'
        write(*,*) '-h                   Print short help'
        write(*,*) '-hall                Print full help'
        write(*,*) '-b <BETA>            Forces the use of beta = <BETA> instead of the'
        write(*,*) '                         value given in the xyz file'
        write(*,*) '--beta <BETA>        Same as above'
        write(*,*) '--data <STRING>      Read bead data from location or the library.'
        write(*,*) '                         See DATA HANDLING below for details.'
        write(*,*) '--pes <PES>          Defines the potential energy surface to use'
        write(*,*) '                         (Only required if LINPOL_DEFAULT_PES'
        write(*,*) '                          not set)'
        write(*,*) '--eps <EPS>          Convergence criterion for L-BFGS search (Ha/bohr)'
        write(*,*) '                         (DEFAULT linpol:      ' , choseneps, ')'
        write(*,*) '                         (DEFAULT single bead: ' , starteps, ')'
        write(*,*) '--gradeps <GRADEPS>  Step size for finite-difference gradient (bohr)'
        write(*,*) '                         Ignored if gradient is exact (e.g. for MB-pol).'
        write(*,*) '                         (DEFAULT: ',gradeps,')'
        write(*,*) '--hesseps <HESSEPS>  Step size for finite-difference second derivative (bohr)'
        write(*,*) '                         (DEFAULT: ',hesseps,')'
        write(*,*) '--writedata <FILE>   Write bead data to <FILE>'
        write(*,*) '                         (Ignored if more than one bead in <INPUT>)'
        write(*,*) '--loose              Use loose linear polymer ends', loosedefault
        write(*,*) '--fixed              Use fixed linear polymer ends', fixeddefault
        write(*,*) '--maxiter <INT>      Maximum number of LBFGS iterations' 
        write(*,*) '                         (DEFAULT: ', lbfgs_maxiter, ')'
        write(*,*) '--noiter             Entirely suppress running of any LBFGS iteration' 
        write(*,*) '--nosplit            Do not calculate Phi (ratio between tunnelling '
        write(*,*) '                          and non-tunnelling determinants) and'
        write(*,*) '                          h (kink path weight)'
        write(*,*) '--nohess             Hessian is not set up or diagonalised,' 
        write(*,*) '                         implies --nosplit' 
        write(*,*) '--debug              Print debug comments'

        if (.not. verbous ) return

        write(*,*) ''
        write(*,*) '====================================================================='
        write(*,*) ''
        write(*,*) 'General description of what the program does'
        write(*,*) '1. read input file'
        write(*,*) '2. optimise instanton using L-BFGS (suppress via --noiter)'
        write(*,*) '3. Calculate Skink and pathway'
        write(*,*) '4. Calculate Hessian, diagonalise it (supress by giving '
        write(*,*) '      switch --nohess)'
        write(*,*) '5. Calculate splitting parameter: logdet, Phi, h (supress by either'
        write(*,*) '      giving switch --nohess or --nosplit)'
        write(*,*) ''
        write(*,*) ''
        write(*,*) 'UNITS AND INPUT FILE CONVENTION'
        write(*,*) 'The xyz file fed into the program should be in Angstrom.  The program'
        write(*,*) 'will then convert it to Bohr for internal use, as Linpol uses atomic'
        write(*,*) 'units (Bohr, Ha, etc.) for all of its calculations.  This convention is'
        write(*,*) 'based on that used by the potentials produced by Joel Bowman and co-workers,'
        write(*,*) 'and they describe the input file requirements as follows:'
        write(*,*) '  "By default, the unit of the Cartesian coordinates is assumed to be'
        write(*,*) '  Angstrom.  Atoms need to be in the order of ascending mass, so "H" always'
        write(*,*) '  goes before "O". All "H" atoms and "O" atoms should be listed in the same'
        write(*,*) '  monomer order.  For example, a water dimer should look like, "H1 H2 H3 H4'
        write(*,*) '  O5 O6", where H1-O5-H2 is one monomer and H3-O6-H4 is the second monomer."'
        write(*,*) ''
        write(*,*) 'Care should be taken when adding new potentials to the program, especially if'
        write(*,*) 'they use different input file conventions and/or different units.'
        write(*,*) 'The pes_interface.f90 source file contains guidance on how to do this.'
        write(*,*) ''
        write(*,*) ''
        write(*,*) 'DATA HANDLING:'
        write(*,*) 'For the proper calculation of Skink and Phi of an instanton the program'
        write(*,*) 'needs information about the local minima this instanton connects.'
        write(*,*) 'A local library with some precalculated data files can be found in the'
        write(*,*) 'folder "data/<PES>/cluster". (where <PES> is the potential used)'
        write(*,*) ''
        write(*,*) 'The file from which this bead data should be read can be selected by'
        write(*,*) 'providing ione of the following  after the flag "--data":'
        write(*,*) '  . the path to a suitable .dat file'
        write(*,*) '  . the relative path to a file within the folder data/<PES>/cluster'
        write(*,*) '        (eg. water/10/PP1.dat)'
        write(*,*) '  . the relative path to a file within a default cluster folder'
        write(*,*) '        (see below)'
        write(*,*) '  . in each case the .dat ending may be omitted'
        write(*,*) ''
        write(*,*) 'If many calculations are done with the same cluster one can save a'
        write(*,*) 'cluster folder to the environment variable LINPOL_DEFAULT_CLUSTER.'
        write(*,*) 'E.g. if this variable is set to "water/10" then "water/10/PP1.dat"'
        write(*,*) 'can be invoked by "PP1.dat" or even just "PP1". If however '
        write(*,*) 'LINPOL_DEFAULT_CLUSTER="water", we would need to provide at least "10/PP1".'
        write(*,*) 'In case <INPUT> has only one bead the --data flag is ignored. '
        write(*,*) 'If <INPUT> has multiple beads and --data is not provided a warning'
        write(*,*) 'will be printed.'
        write(*,*) ''
        write(*,*) 'Suitable data files may be created by --writedata <FILE>.'
        write(*,*) 'This only works if <INPUT> only has a single bead and neither "--nohess", '
        write(*,*) 'nor "--noiter" were given. Otherwise the flag is ignored.'
        write(*,*) ''
        write(*,*) ''
        write(*,*) 'ENVIRONMENT VARIABLES:'
        write(*,*) 'LINPOL: has to be set to the folder containing the linpolexe'
        write(*,*) 'LINPOL_DEFAULT_CLUSTER: may be set to a folder within data/<PES>/cluster'
        write(*,*) 'LINPOL_DEFAULT_PES: may be set to a default PES'
        write(*,*) 'You may also need to set other environment variables for the potentials,'
        write(*,*) 'but it should become apparent which are required when you try running the'
        write(*,*) 'program for the first time.  For example, HBB2pol requires an environment'
        write(*,*) 'variable pointing to its data files.'
        write(*,*) ''
        write(*,*) ''
        write(*,*) 'GENERAL NOTES:'
        write(*,*) 'The format for <EPS>, <GRADEPS> and <HESSEPS> is 0.0D0.'
        write(*,*) 'Decreasing <GRADEPS> might have only a small influence on the tightness'
        write(*,*) '    of convergence, more on the speed. Too small values will, however,'
        write(*,*) '    be likely to cause numerical problems. Some experiments have shown'
        write(*,*) '    that best value for <GRADEPS> is in the region of 1E-4 to 1E-3.'
        write(*,*) 'The optimal value of <HESSEPS> is usually around 2.0D-2, which is'
        write(*,*) '    sufficiently small to produce a good approximation to the second'
        write(*,*) '    derivative, but sufficiently large for the calculation not to suffer'
        write(*,*) '    from small fluctuations in the potential between closely-spaced points.'
        write(*,*) 'When --eps is set both the values for single bead and the linear polymer'
        write(*,*) '    convergence are altered. If --eps is not touched convergence is '
        write(*,*) '    to ',starteps,' in single bead'
        write(*,*) '    and to ', choseneps,' in linpol mode.'
        write(*,*) 'The last argument of --loose, --fix dominates.'
        write(*,*) 'Fixed ends requires files "minA" and "minB" giving the positions of the minima.'
        write(*,*) 'These files should contain only the coordinates for the atoms; they do not need'
        write(*,*) 'to be in full xyz format.  Indeed, xyz files will not be read successfully.'
end subroutine

function valid_file(f)
!test if a file f is actually a valid, ie existing and readable file
CHARACTER(len=*),INTENT(IN) :: f
logical :: valid_file


CHARACTER(80) :: msg                    ! Status message
INTEGER :: stat                         ! IO status message

OPEN (UNIT=999, FILE=f, STATUS='OLD', ACTION='READ', IOSTAT=stat, IOMSG=msg )
valid_file = (stat .eq. 0)
CLOSE (UNIT=999)
end function

subroutine read_double(posn,switch,tovar)
!test for consistency and read a double from position <posn> following the
!switch <switch>, saving the result to <tovar>
        integer,intent(in) :: posn
        character(len=*), intent(in) :: switch
        character(len=512) :: tmp
        real(kind=dp), intent(out) :: tovar
        INTEGER :: stat                         ! IO status message

        if (posn .ge. argc) then
                if (myrank .eq. root) write(*,*) 'The switch "',trim(switch),'" needs a double argument following'
                stop 1
        end if
        call getarg(posn,tmp)
        READ(tmp,*,IOSTAT=stat) tovar
        if (stat .ne. 0) then
                if (myrank .eq. root) write(*,*) 'Please provide a valid floating point value after "',trim(switch),'".'
                stop 1
        end if
end subroutine

subroutine read_integer(posn,switch,tovar)
!test for consistency and read an int from position <posn> following the
!switch <switch>, saving the result to <tovar>
        integer,intent(in) :: posn
        character(len=*), intent(in) :: switch
        character(len=512) :: tmp
        integer, intent(out) :: tovar
        INTEGER :: stat                         ! IO status message

        if (posn .ge. argc) then
                if (myrank .eq. root) write(*,*) 'The switch "',trim(switch),'" needs an integer argument following'
                stop 1
        end if
        call getarg(posn,tmp)
        READ(tmp,*,IOSTAT=stat) tovar
        if (stat .ne. 0) then
                if (myrank .eq. root) write(*,*) 'Please provide a valid integer value after "',trim(switch),'".'
                stop 1
        end if
end subroutine

subroutine read_character(posn,switch,tovar)
!test for consistency and read a char from position <posn> following the
!switch <switch>, saving the result to <tovar>
        integer,intent(in) :: posn
        character(len=*), intent(in) :: switch
        character(len=512),intent(out) :: tovar
        INTEGER :: stat                         ! IO status message

        if (posn .ge. argc) then
                if (myrank .eq. root) write(*,*) 'The switch "',trim(switch),'" needs a string argument following'
                stop 1
        end if
        call getarg(posn,tovar)
end subroutine

subroutine read_pes_type(posn,switch,tovar)
!test for consistency and read a char from position <posn> following the
!switch <switch>, saving the result to <tovar>
        integer,intent(in) :: posn
        character(len=*), intent(in) :: switch
        character(len=512) :: tmp
        integer,intent(out) :: tovar
        INTEGER :: stat                         ! IO status message

        if (posn .ge. argc) then
                if (myrank .eq. root) write(*,*) 'The switch "',trim(switch),'" needs a <PES> argument following'
                stop 1
        end if
        call getarg(posn,tmp)
        call pes_type_from_string(tmp,tovar)

        if (tovar .eq. PES_NONE) then
                if (myrank .eq. root) write(*,*) 'Please provide a valid <PES> name after "',trim(switch),'".'
                write(*,*) "The following potentials are implemented:"
                
                tovar=1
                do while (.true.)
                       call pes_string_from_type(tovar,tmp)

                       if (trim(tmp) .ne. "NONE") then
                               write(*,*) "    - ",trim(tmp)
                       else
                               exit
                       end if
                       tovar = tovar+1
                end do
                stop 1
        end if
end subroutine

subroutine parse_args()
        use mpi
        integer :: i,stat
        character(len=512) :: tmp         !Char of undetermined length

        ! fill modular vars:
        argc = iargc()
        CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)
        filename=""     !initialise as empty

        !-----------------------

        if ( argc .lt. 1 ) then
                if (myrank .eq. root) write(*,*) 'Please provide at least one argument (the filename)'
                stop 1
        endif  

        i=1
        do
                CALL getarg ( i, tmp )
                select case(tmp)
                        case ("-h")
                                call print_help(.false.)
                                stop
                        case ("-hall")
                                call print_help(.true.)
                                stop
                        case ("-b")
                                i=i+1
                                call read_double(i,"-b",masterbeta)
                        case ("--beta")
                                i=i+1
                                call read_double(i,"--beta",masterbeta)
                        case ("--writedata")
                                i=i+1
                                call read_character(i,"--writedata",dataoutputfile)
                        case ("--data")
                                i=i+1
                                call read_character(i,"--data",datainput)
                                !NOTE: this file is checked for validity only if
                                !neccessary in linpolexe.f90
                        case ("--pes")
                                i=i+1
                                call read_pes_type(i,"--pes",pes_type)
                        case ("--loose")
                                endchoice=2
                        case ("--fixed")
                                endchoice=1
                        case ("--debug")
                                debug = .true.
                        case ("--maxiter")
                                i=i+1
                                call read_integer(i,"--maxiter",lbfgs_maxiter)
                                if (lbfgs_maxiter .le. 0) then
                                        do_iter = .false.
                                end if
                        case ("--eps")
                                i=i+1
                                call read_double(i,"--eps",choseneps)
                                starteps=choseneps      !copy value to starteps (used for single bead mode)
                        case ("--gradeps")
                                i=i+1
                                call read_double(i,"--gradeps",gradeps)
                        case ("--hesseps")
                                i=i+1
                                call read_double(i,"--hesseps",hesseps)
                        case ("--noiter")
                                do_iter = .false.
                        case ("--nosplit")
                                do_split = .false.
                        case ("--nohess")
                                do_split = .false.
                                do_hess = .false.
                        case default
                                if ( i .ge. argc ) then
                                        if (valid_file(tmp)) then
                                                filename = tmp
                                        else
                                                if (myrank .eq. root) &
                                                        write(*,*) "Please provide a valid xyz input file, not ",trim(tmp)
                                                stop 1
                                        end if
                                else
                                        if (myrank .eq. root) &
                                                write (*,*) "Unknown switch: ",trim(tmp)
                                        stop 1
                                end if
                end select

                if (i .lt. argc) then
                        i=i+1
                else
                        exit ! exit loop
                endif
        enddo

        !do we actually have a valid file name?
        if (filename .eq. '') then
                if (myrank .eq. root) write(*,*) "Please provide a valid xyz input file, not ",trim(filename)
                stop 1
        end if
        if (pes_type .eq. PES_NONE) then
                CALL GET_ENVIRONMENT_VARIABLE("LINPOL_DEFAULT_PES", tmp)
                call pes_type_from_string(tmp,pes_type)

                if (pes_type .eq. PES_NONE) then
                        if (myrank .eq. root) then
                                write(*,*) "Please provide the PES to be used via --pes <PES>"
                                write(*,*) "or by setting LINPOL_DEFAULT_PES"
                                write(*,*) "The following potentials are implemented:"
                                
                                i=1
                                do while (.true.)
                                       call pes_string_from_type(i,tmp)

                                       if (trim(tmp) .ne. "NONE") then
                                               write(*,*) "    - ",trim(tmp)
                                       else
                                               exit
                                       end if
                                       i=i+1
                                end do
                        end if
                        stop 1
                end if
        end if
end subroutine

end module
