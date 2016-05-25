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
PROGRAM linpolexe
USE funcs						! Miscellaneous functions module
USE prep						! Linear polymer preparation module
USE pots						! Potential calculation module
USE grads						! Gradient calculation module
USE hess						! Hessian calculation module
USE lbfgs_module					! Marko's L-BFGS module
USE parameters						! Module containing parameters to be pulled into all other modules
USE clusterdatamod                                      ! contains the cluster data strucuture and functions
USE min							! Module containing subroutines for locating minima
USE mpi							! Load MPI for parallelization
use pes_interface
use opt_interface
use args
use sort
IMPLICIT NONE
!------------

! =======================================================================
! PURPOSE OF THIS PROGRAM
! =======================================================================
! This program was written to calculate the linear-polymer potential for a
! system of beads described by an xyz file.
!
! Record of revisions:
!
!  Date           Person         Change implemented
!  ----           ------         ------------------
!  14 Sep 2011    A. Reid        Original code developed
!  Oct 2011       A. Reid        Ongoing development
!  27 Oct 2011    A. Reid        Splitting of code into subroutines
!  1 Nov 2011     A. Reid        Adding gradient subroutines
!  4 Nov 2011     A. Reid        Adding start/end and beadize subroutines
!  23 Nov 2011    A. Reid        Integration with MTC's minimization code
!  29 Nov 2011    A. Reid        Conversion to use variables from modules
!  9 Jan 2012     A. Reid        Interpolation deferred to Python
!  12 Jan 2012    A. Reid        Finished debugging the two modes
!  14 Feb 2012    A. Reid        Changed Hessian to banded storage
!  8 May 2012     A. Reid        Removed old commented-out code
!  8 May 2012     A. Reid        Removed interpolation code (old 'case 2')
!  8 May 2012     A. Reid        Format descriptors tidied up and IF comments removed
!  15 May 2012    A. Reid        Code reorganised into different cases (1 and 2)
!  11 Jun 2012    A. Reid        Parallelized with MPI
!  27 Jun 2012    A. Reid        Fixed single bead case
!  20 Nov 2012    M. Herbst      Implemented command line args
!  24 Nov 2012    M. Herbst      Parallelised hessian calculation routines
!  23 Dec 2012    M. Herbst      Split linpolexe into subroutines
!
! =======================================================================

! =======================================================================
! DECLARATION OF VARIABLES
! =======================================================================
INTEGER :: error, nprocs, namelen, myrank		! For MPI
!CHARACTER ( LEN = MPI_MAX_PROCESSOR_NAME ) :: procname

double precision :: startt, mint, hesst, hessDiagT         ! for duration calc
double precision,DIMENSION(:),ALLOCATABLE :: evalues       ! to hold evalues
double precision :: skink  ! the value for skink
!CHARACTER ( LEN = 128 ) :: hessianoutput = 'hess.out'	! File for Hessian output (used for debugging purposes)

! =======================================================================
! INITIALIZE MPI AND RECORD MPI START TIME
! =======================================================================
CALL MPI_Init(error)
CALL MPI_Comm_size ( MPI_COMM_WORLD, nprocs, error)
CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)
!CALL MPI_Get_processor_name (procname, namelen, error)
startt = MPI_Wtime()
mint = 0
hesst = 0
hessDiagT = 0

! =======================================================================
! READ COMMAND LINE ARGS
! =======================================================================
call parse_args()

! =======================================================================
! SET UP FILENAMES USING DATE AND TIME AND START LOG FILE
! =======================================================================
!filenames are stored in parameter.f90
call do_setup_output_files()

! =======================================================================
! LOG INTRODUCTORY TEXT
! =======================================================================
call do_start_log()     !NOTE: unit 99 is log file
call do_introduction()

! =======================================================================
! INITIAL READ OF FILE TO ESTABLISH SIZE
! =======================================================================
CALL initialends ( filename )
n = ninit ! Set no. of beads equal to value read in from xyz file

! =======================================================================
! CHECK FOR SANITY OF MPI ENVIRONMENT
! =======================================================================
call do_check_mpi_sanity(ninit)

! =======================================================================
! INITIALISE THE POTENTIALS
! =======================================================================
call do_init_modules()
CALL MPI_Barrier(MPI_COMM_WORLD, error) !synchronises all processes


! =======================================================================
! SPLIT INTO DIFFERENT MODES
! =======================================================================
if (ninit == 1) then
        !single bead mode
        call do_read_input_single()

        if (do_iter) then
                call do_minimisation_single()
                !first order property calculation:
                call do_property_calculation_single(.false.)    
                !no need to calculate the potential again => false
        else
                call do_property_calculation_single(.true.)
        end if

        if (do_hess) then
                call do_hessian_single(evalues)

                !save data to clusterfile if so whished:
                if (trim(dataoutputfile) .ne. '') then
                        if (.not. writeClusterData(dataoutputfile,size(evalues),latest_gnorm,getdatetime(),singlepot,evalues)) then
                                write (*,*) "ClusterData file could not be written !"
                        end if
                end if
        end if
else
        !multiple bead mode
        call do_read_input_multiple()

        if (do_iter) then
                call do_minimisation_multiple()
                !first order property calculation:
                call do_property_calculation_multiple(.false.,skink)
                !no need to calculate the potential again => false
        else
                call do_property_calculation_multiple(.true.,skink)
        end if

        if (do_hess) then
                call do_hessian_multiple(evalues)

                if (do_split) then
                        call do_splitting_multiple(evalues,skink)
                end if
        end if
end if
call do_end(.true.) !end the program for a successful run

!#############################################################################################################################
!#############################################################################################################################

contains
        subroutine do_setup_output_files()
                ! The following filenames are common to both modes (single and multiple) of the program
                initoutputfile = getdatetime()//'-unoptimized.xyz'
                outputfile = getdatetime()//'-optimized.xyz'
                logfile = getdatetime()//'-log.txt'
                pathwayfile = getdatetime()//'-pathway.txt'
                savefile = getdatetime()//'-saved.xyz'
                lowestgradfile = getdatetime()//'-lowestgrad.xyz'
                lowestfuncfile = getdatetime()//'-lowestfunc.xyz'
        end subroutine

!=============================================================================================================================

        subroutine do_start_log()
                ! Create a single log file in the root process
                IF (myrank == root) THEN
                        OPEN(UNIT=99,FILE=logfile,STATUS='new',ACTION='write')
                END IF
        end subroutine

!=============================================================================================================================

        subroutine do_introduction()
                !log introductory stuff
                if (myrank .ne. root) return

                ! =======================================================================
                ! INTRODUCTORY TEXT
                ! =======================================================================
                WRITE(99,'(A)')
                WRITE(99,'(A)') '=======================================================================' 
                WRITE(99,'(A)') 'LINEAR-POLYMER INSTANTON CODE'
                WRITE(99,'(A)') '======================================================================='
                WRITE(*,'(A)')
                WRITE(*,'(A)') '=======================================================================' 
                WRITE(*,'(A)') 'LINEAR-POLYMER INSTANTON CODE'
                WRITE(*,'(A)') '======================================================================='

!		WRITE(99,'(A)') ''
!		WRITE(99,'(A)') '##################################'
!		WRITE(99,'(A)') 'USING FIVE-POINT STENCIL GRADIENT!'
!		WRITE(99,'(A)') '##################################'
!		WRITE(99,'(A)') ''
!		WRITE(*,'(A)') ''
!		WRITE(*,'(A)') '##################################'
!		WRITE(*,'(A)') 'USING FIVE-POINT STENCIL GRADIENT!'
!		WRITE(*,'(A)') '##################################'
!		WRITE(*,'(A)') ''

                ! =======================================================================
                ! REPORT FILENAMES
                ! =======================================================================

                WRITE(99,'(A)') 'XYZ file chosen by user to provide unoptimized geometry = ', filename
                WRITE(99,'(A)') 'XYZ file containing copy of unoptimized geometry = ', initoutputfile
                WRITE(99,'(A)') 'XYZ file containing optimized geometry (post L-BFGS) = ', outputfile
                WRITE(99,'(A)') 'Text file containing log of program execution = ', logfile
                WRITE(*,'(A)') 'Text file containing log of program execution = ', logfile
        end subroutine

!=============================================================================================================================

        subroutine do_check_mpi_sanity(ninit)
                integer, intent(in) :: ninit !the number of beads in input file
                INTEGER :: nbeads            ! Number of beads allocated to each process

                if (myrank .ne. root) then
                        !do a quick check and exit if wrong:
                        if (mod(ninit,nprocs) .ne. 0) then
                                call do_end(.false.) ! end without successful calculation
                        end if
                        return
                end if

                !we are root process:
                WRITE(*,'(A)') '------------------------------------------------------------------'
                WRITE(*,'(A)') 'PARALLELIZATION: ALLOCATION OF BEADS ACROSS PROCESSES'
                WRITE(*,'(A)') '------------------------------------------------------------------'
                WRITE(99,'(A)') '------------------------------------------------------------------'
                WRITE(99,'(A)') 'PARALLELIZATION: ALLOCATION OF BEADS ACROSS PROCESSES'
                WRITE(99,'(A)') '------------------------------------------------------------------'

                IF (nprocs .GT. ninit) THEN ! A. More processes than beads
                        WRITE(*,'(A,I5,A,I5,A)') 'There are ', ninit, ' bead(s) and ', nprocs, ' processes;'
                        WRITE(*,'(A)') 'ERROR: No. of processes is greater than no. of beads!'
                        WRITE(*,'(A)') 'The program will now exit.'
                        WRITE(99,'(A,I5,A,I5,A)') 'There are ', ninit, ' bead(s) and ', nprocs, ' processes;'
                        WRITE(99,'(A)') 'ERROR: No. of processes is greater than no. of beads!'
                        WRITE(99,'(A)') 'The program will now exit.'
                        call do_end(.false.) ! end without successful calculation
                        return
                end if

                if (mod(ninit,nprocs) .ne. 0) then ! C. Number of processes is not integer factor of number of beads
                        WRITE(*,'(A,I5,A,I5,A)') 'There are ', ninit, ' beads and ', nprocs, ' processes;'
                        WRITE(*,'(A)') 'Hence, the beads cannot be allocated evenly across the processes.'
                        WRITE(*,'(A)') 'ERROR: No. of processes is not an integer factor of the no. of beads!'
                        WRITE(*,'(A)') 'The program will now exit.'
                        WRITE(99,'(A,I5,A,I5,A)') 'There are ', ninit, ' beads and ', nprocs, ' processes;'
                        WRITE(99,'(A)') 'Hence, the beads cannot be allocated evenly across the processes.'
                        WRITE(99,'(A)') 'ERROR: No. of processes is not an integer factor of the no. of beads!'
                        WRITE(99,'(A)') 'The program will now exit.'
                        call do_end(.false.) ! end without successful calculation
                        return
                end if

                nbeads = ninit/nprocs
                WRITE(*,'(A,I5,A,I5,A)') 'There are ', ninit, ' bead(s) and ', nprocs, ' process(es);'
                WRITE(*,'(A,I5,A)') 'Hence, ', nbeads, ' bead(s) will be allocated to each process.'
                WRITE(99,'(A,I5,A,I5,A)') 'There are ', ninit, ' bead(s) and ', nprocs, ' process(es);'
                WRITE(99,'(A,I5,A)') 'Hence, ', nbeads, ' bead(s) will be allocated to each process.'

                WRITE(*,'(A)') '------------------------------------------------------------------'
                WRITE(99,'(A)') '------------------------------------------------------------------'
        end subroutine

!=============================================================================================================================
        
        subroutine do_init_modules()
                character(len=8) :: tmp

                call pes_obj%init(pes_type, p/3)
                call opt_obj%init(opt_type)

                if (myrank .eq. root) then
                        call pes_obj%pes_name(tmp)
                        write(*,*)
                        write(99,*)
                        write(*,'(2A)') "Using PES: ", tmp
                        write(99,'(2A)') "Using PES: ", tmp

                        call opt_obj%opt_name(tmp)
                        write(*,'(2A)') "Using OPT method: ", tmp
                        write(99,'(2A)') "Using OPT method: ", tmp
                        write(*,*)
                        write(99,*)
                end if
        end subroutine

!=============================================================================================================================

        subroutine do_read_input_single()
                ! =======================================================================
                ! SET UP MASSES ARRAY, AND ATOM TAGS ARRAY
                ! =======================================================================
                CALL singlebead ( filename )
        end subroutine

!=============================================================================================================================

        subroutine do_read_input_multiple()
                logical :: success
                ! =======================================================================
                ! SET UP STARTING AND ENDING ARRAYS, MASSES ARRAY, AND ATOM TAGS ARRAY
                ! =======================================================================

                !read cluster data:
                if (myrank .eq. root) then
                        WRITE(*,'(A)') ''
                        WRITE(*,'(A)') 'READING CLUSTER DATA'
                        WRITE(99,'(A)') ''
                        WRITE(99,'(A)') 'READING CLUSTER DATA'
                end if

                ! error checking on datainput can only be done now that we know
                ! that we deal with a multiple bead case!
                if (trim(datainput) .eq. '') then
                        if (myrank .eq. root) then
                                write(*,*) "WARNING: no cluster data file given."
                                write(*,*) "Potential is not offset: Skink values not correctly shifted"
                                write(*,*) "Calculation of phi will be suppressed"
                                write(99,*) "WARNING: no cluster data file given."
                                write(99,*) "Potential is not offset: Skink values not correctly shifted"
                                write(99,*) "Calculation of phi will be suppressed"
                        end if
                        do_split = .false.
                        cdata%potential = 0.0
                else if (is_valid_cluster_data_string(datainput,pes_obj)) then
                        call readClusterData(cdata,datainput,pes_obj,success)
                        if (success) then !read data to container => true all went well, false: things went wrong
                                !all went well

                                if (myrank .eq. root) then
                                        WRITE(99,'(A)') '-----------------------------------------------------------------------'
                                        write(99,'(A,a)') 'Path of file read: ', trim(cdata%fullPath)
                                        WRITE(99,'(A)') 'Some of the data read:'
                                        write(99,'(6x,A,i20)')    'degrees of freedom = ', cdata%dof
                                        write(99,'(6x,A,E20.13)') 'gnorm              = ', cdata%gnorm
                                        write(99,'(6x,a,a)')      'date-time          = ', cdata%dateTime
                                        write(99,'(6x,a,a)')      'comments           = ', trim(cdata%comments)
                                        write(99,'(6x,A,E20.13)') 'potential          = ', cdata%potential
                                        if (cdata%global_minimum) then
                                                write(99,'(6x,a)') 'This is the GLOBAL MINIMUM'
                                        end if
                                        write(99,*)


                                        WRITE(*,'(A)') '-----------------------------------------------------------------------'
                                        WRITE(*,'(A)') 'Some of the data read:'
                                        write(*,'(6x,A,i20)')    'degrees of freedom = ', cdata%dof
                                        write(*,'(6x,A,E20.13)') 'gnorm              = ', cdata%gnorm
                                        write(*,'(6x,A,E20.13)') 'potential          = ', cdata%potential
                                        write(*,*)
                                end if
                        else
                                if (myrank .eq. root) then
                                        WRITE(99,'(A)') 'Error occurred when reading cluster data file.'
                                        WRITE(99,*) TRIM(get_cluster_data_full_path(datainput,pes_obj,success))
                                        write(99,*) "Check that the runtime variables are set up properly."
                                        WRITE(*,'(A)') 'Error occurred when reading cluster data file.'
                                        WRITE(*,*) trim(get_cluster_data_full_path(datainput,pes_obj,success))
                                        write(*,*) "Check that the runtime variables are set up properly."
                                end if
                                stop 1
                        end if
                else
                        if (myrank .eq. root) then
                                WRITE(99,'(A)') 'An invalid string was passed to program after --data.'
                                WRITE(99,'(A)') 'Could not read cluster data.'
                                WRITE(*,*) 'An invalid string was passed to program after --data.'
                                WRITE(*,*) 'Could not read cluster data.'
                        end if
                        stop 1
                end if

                CALL array ( filename )

                if (endchoice==1) then
                        allocate(xA(d,p), xB(d,p))
                        open(21,file='minA')
                        read(21,*) xA
                        close(21)
                        open(21,file='minB')
                        read(21,*) xB
                        close(21)
                        !convert from Angstroms to Bohr
                        xA = xA / auang
                        xB = xB / auang
                end if

                CALL data ( )    ! REPORTING DATA ON 3-D ARRAY AND CALCULATING BETAN
                
        end subroutine

!=============================================================================================================================

        subroutine do_minimisation_single()
                if (myrank .ne. root) return      !don't do a thing if we are not root

                WRITE(99,'(A)')
                WRITE(99,'(A)') '-----------------------------------------------------------------------'
                WRITE(99,'(A)') 'CALCULATING OPTIMUM GEOMETRY FOR SINGLE BEAD'
                WRITE(99,'(A)') '-----------------------------------------------------------------------'
                WRITE(*,'(A)') '-----------------------------------------------------------------------'
                WRITE(*,'(A)') 'CALCULATING OPTIMUM GEOMETRY FOR SINGLE BEAD'
                WRITE(99,'(A)')

                ! =======================================================================
                ! PREPARE TO FIND MINIMUM ENERGY GEOMETRY
                ! =======================================================================
                CALL xyzout ( xsingle , n , d , p, beta, atomtags, initoutputfile )
                WRITE(99,'(A)') 'Unoptimized bead geometry written to initoutputfile'

                ! =======================================================================
                ! CALL MINIMIZATION ROUTINE AND THEN REPORT RESULTS
                ! =======================================================================
                WRITE(99,'(A)') '-----------------------------------------------------------------------'
                WRITE(99,'(A)') "Running the L-BFGS code to find cluster geometry"
                WRITE(99,'(A)') '-----------------------------------------------------------------------'
                WRITE(*,'(A)') "Running the L-BFGS code to find cluster geometry"
                WRITE(*,'(A)') '-----------------------------------------------------------------------'

                mint = MPI_Wtime()
                CALL minsingle ( outputfile )
                mint = MPI_Wtime() - mint

                WRITE(99,'(A)') '-----------------------------------------------------------------------'
                WRITE(99,'(A)') "Minimization finished (either converged, or maximum iterations reached)"
                WRITE(99,'(A)') '-----------------------------------------------------------------------'
                WRITE(*,'(A)') "Minimization finished (either converged, or maximum iterations reached)"
                WRITE(*,'(A)') '-----------------------------------------------------------------------'
        end subroutine

!=============================================================================================================================

        subroutine do_minimisation_multiple() 

                ! =======================================================================
                ! PRINT OUT INFORMATION FOR USER AND WRITE UNOPTIMIZED FILE
                ! =======================================================================
                IF (myrank == root) THEN			! To avoid duplication in non-root processes
        
                        WRITE(99,'(A)')
                        WRITE(99,'(A)') '-----------------------------------------------------------------------'
                        WRITE(99,'(A)') 'CALCULATING OPTIMUM GEOMETRY FOR LINEAR POLYMER'
                        WRITE(99,'(A)') '-----------------------------------------------------------------------'
                        WRITE(*,'(A)') '-----------------------------------------------------------------------'
                        WRITE(*,'(A)') 'CALCULATING OPTIMUM GEOMETRY FOR LINEAR POLYMER'

                        CALL xyzout ( x , n , d , p, beta, atomtags, initoutputfile )
                        WRITE(99,'(A)') 'Unoptimized x written to initoutputfile'

                        WRITE(99,'(A)') '-----------------------------------------------------------------------'
                        WRITE(99,'(A)') "Running the L-BFGS minimization code to find instanton"
                        WRITE(99,'(A)') '-----------------------------------------------------------------------'
                        WRITE(*,'(A)') "Running the L-BFGS minimization code to find instanton"
                        WRITE(*,'(A)') '-----------------------------------------------------------------------'

                end if

                ! =======================================================================
                ! CALL MINIMIZATION ROUTINE AND THEN REPORT RESULTS
                ! =======================================================================

                mint = MPI_Wtime()
                CALL minlinpol ( outputfile )
                mint = MPI_Wtime() -mint

                IF (myrank == root) THEN			! To avoid duplication in non-root processes
                        WRITE(99,'(A)') '-----------------------------------------------------------------------'
                        WRITE(99,'(A)') "Minimization finished (either converged, or maximum iterations reached)"
                        WRITE(99,'(A)') '-----------------------------------------------------------------------'
                        WRITE(*,'(A)') '-----------------------------------------------------------------------'
                        WRITE(*,'(A)') "Minimization finished (either converged, or maximum iterations reached)"
                        WRITE(*,'(A)') '-----------------------------------------------------------------------'
                END IF					! To avoid duplication in non-root processes
        end subroutine

!=============================================================================================================================

        subroutine do_property_calculation_single(calc_potential_gradient)
                ! calculate first order properties here
                ! ie the ones involving the potential and the gradient
                ! but not hessian or higher derivatives

                logical, intent(in) :: calc_potential_gradient
                !      true: the global singlepot and latest_gradient vars have not been filled and
                !            hence we need to do that
                !      false: it is filled and we can use var straight away


                WRITE(99,'(A)') '######################################################################'
                WRITE(99,'(A)') "#---                  PROPERTY CALCULATION                        ---#"
                WRITE(99,'(A)') '######################################################################'
                WRITE(*,'(A)') '######################################################################'
                WRITE(*,'(A)') "#---                  PROPERTY CALCULATION                        ---#"
                WRITE(*,'(A)') '######################################################################'

                WRITE(99,'(A,ES10.3,A)') ' GRADEPS for finite-difference gradient set to ', gradeps, ' bohr'
                WRITE(*,'(A,ES10.3,A)') ' GRADEPS for finite-difference gradient set to ', gradeps, ' bohr'

                if (calc_potential_gradient) then
                        !##########################################
                        !#-- transform beads to pes coordinates --#
                        !##########################################
                        call pes_obj%atomvec_to_pes(masses)
                        call pes_obj%vec_to_pes(xsingle)
                        !call pes_obj%vec_to_pes(atomtags)
                        transformed_to_pot = .true.
                        !##########################################

                        call pes_obj%pot(xsingle,singlepot)
                        latest_gnorm = onebeadgnorm(reshape(xsingle, (/ d*p /) ))

                        !############################################
                        !#-- transform beads from pes coordinates --#
                        !############################################
                        call pes_obj%vec_from_pes(xsingle)
                        call pes_obj%atomvec_from_pes(masses)
                        !call pes_obj%vec_from_pes(atomtags)
                        transformed_to_pot = .false.
                        !############################################
                end if

		! These values are always in Ha and Ha/bohr regardless of the potential (see pes_interface.f90)
                WRITE(99,*) 'Potential: ', singlepot, ' hartree'
                WRITE(*,*) 'Potential: ', singlepot, ' hartree'
                WRITE(99,*) 'GNORM: ', latest_gnorm, ' hartree/bohr'
                WRITE(*,*) 'GNORM: ', latest_gnorm, ' hartree/bohr'
        end subroutine

!=============================================================================================================================
        
        subroutine do_property_calculation_multiple(calc_potential,skink)
                ! calculate first order properties here
                ! ie the ones involving the potential and the gradient
                ! but not hessian or higher derivatives

                logical, intent(in) :: calc_potential   
                !      true: the global totalpot var has not been filled and
                !            hence we need to do that
                !      false: it is filled and we can use var straigt away
                double precision,intent(out) :: skink					! Skink, as calculated with the selected 'zeropot' value (see parameters.f90)

                IF (myrank == root) THEN			! To avoid duplication in non-root processes
                        WRITE(99,'(A)') ''
                        WRITE(99,'(A)') '######################################################################'
                        WRITE(99,'(A)') "#---                  PROPERTY CALCULATION                        ---#"
                        WRITE(99,'(A)') '######################################################################'
                        WRITE(*,'(A)') ''
                        WRITE(*,'(A)') '######################################################################'
                        WRITE(*,'(A)') "#---                  PROPERTY CALCULATION                        ---#"
                        WRITE(*,'(A)') '######################################################################'
                END IF					! To avoid duplication in non-root processes

                if (calc_potential) then
                        !##########################################
                        !#-- transform beads to pes coordinates --#
                        !##########################################
                        do i=1,n
                                call pes_obj%vec_to_pes(x(:,:,i))
                        end do
                        call pes_obj%atomvec_to_pes(masses)
                        !call pes_obj%vec_to_pes(atomtags) -- should ideally do this, too 
                        ! (Omitted right now since not needed; should be implemented in linpole object, however)
                        if (endchoice == 1) then !fixed ends
                                call pes_obj%vec_to_pes(xA)
                                call pes_obj%vec_to_pes(xB)
                        end if
                        transformed_to_pot = .true.
                        !############################################

                        call totpotx(reshape(x, (/ d*p*n /) ),totalpot)

                        !############################################
                        !#-- transform beads from pes coordinates --#
                        !############################################
                        do i=1,n
                                call pes_obj%vec_from_pes(x(:,:,i))
                        end do
                        call pes_obj%atomvec_from_pes(masses)
                        !call pes_obj%vec_from_pes(atomtags)
                        if (endchoice == 1) then !fixed ends
                                call pes_obj%vec_from_pes(xA)
                                call pes_obj%vec_from_pes(xB)
                        end if
                        transformed_to_pot = .false.
                        !############################################
                end if

                skink = betan * hbar * totalpot
              ! skink = beta/n * hbar * totalpot
              !       = (1/kT)/n * hbar * totalpot
              !       = (1/kTn) * hbar * totalpot
     ! units of skink = (1/k(hartree/k)) * hbar * hartree
     !                = (1/hartree) * hbar * hartree
     !                = hbar
                ! add many more

                if (myrank == root) THEN
                        WRITE(99,*) 'Linear-polymer potential: ', totalpot, ' hartree'
                        WRITE(*,*) 'Linear-polymer potential: ', totalpot, ' hartree'
                        WRITE(99,*) 'Skink: ', skink, ' hbar'
                        WRITE(*,*) 'Skink: ', skink, ' hbar'
                        CALL pathway(x,pathwayfile)
                        WRITE(99,*) 'Instanton pathway saved to:'
                        WRITE(99,*) pathwayfile
                        WRITE(*,*) 'Instanton pathway saved to:'
                        WRITE(*,*) pathwayfile
                END IF
        end subroutine

!=============================================================================================================================

        subroutine do_hessian_single(evalues) !setup and diagonalise hessian matrix for single bead 
                double precision,DIMENSION(:,:),ALLOCATABLE :: hessian
                double precision,DIMENSION(:,:),ALLOCATABLE :: evectors	! Computed eigenvectors of matrix
                double precision,DIMENSION(:),ALLOCATABLE,intent(out) :: evalues	! Computed eigenvalues of matrix

                ! =======================================================================
                ! USE MINIMIZED GEOMETRY TO CALCULATE HESSIAN
                ! =======================================================================

                if (myrank .ne. root) return      !don't do a thing if we are not root

                WRITE(99,'(A)') '-----------------------------------------------------------------------'
                WRITE(99,'(A)') 'Calculating Hessian matrix of second derivatives'
                WRITE(99,'(A)') '-----------------------------------------------------------------------'
                WRITE(*,'(A)') '-----------------------------------------------------------------------'
                WRITE(*,'(A)') 'Calculating Hessian matrix of second derivatives'
                WRITE(*,'(A)') '-----------------------------------------------------------------------'
                WRITE(*,'(A,ES10.3,A)') 'HESSEPS for Hessian calculation set to ', hesseps, ' bohr'
                WRITE(99,'(A,ES10.3,A)') 'HESSEPS for Hessian calculation set to ', hesseps, ' bohr'

                hesst = MPI_Wtime()
                ALLOCATE( hessian((d*p),(d*p)) )	! Create the blank array for the Hessian matrix of second derivatives
                CALL hesscalcsingle ( reshape(xsingle, (/ d*p /) ), hessian )
                hesst = MPI_Wtime() -hesst

                ! =======================================================================
                ! DIAGONALIZE HESSIAN
                ! =======================================================================

                WRITE(99,'(A)') '-----------------------------------------------------------------------'
                WRITE(99,'(A)') 'Diagonalizing Hessian matrix of second derivatives'
                WRITE(99,'(A)') '-----------------------------------------------------------------------'
                WRITE(*,'(A)') '-----------------------------------------------------------------------'
                WRITE(*,'(A)') 'Diagonalizing Hessian matrix of second derivatives'
                WRITE(*,'(A)') '-----------------------------------------------------------------------'

                hessDiagT = MPI_Wtime()
                ALLOCATE( evectors((d*p),(d*p)) )	! Create the blank array for the eigenvectors of the Hessian
                ALLOCATE( evalues(d*p) )		! Create the blank array for the eigenvalues of the Hessian
                CALL diagsymm ( hessian, evalues, evectors )
                hessDiagT = MPI_Wtime() - hessDiagT

                WRITE(99,'(A)')
                WRITE(99,*) 'evalues = '
                WRITE(99,*) evalues
                WRITE(99,'(A)')
                WRITE(*,'(A)')
                WRITE(*,*) 'evalues(:12) = '
                WRITE(*,*) evalues(:12)
                WRITE(*,'(A)')

                ! Vibrational ZPE for harmonic oscillator = (hbar*omega)/2
                ! However, the lowest six eigenvalues have a tendency to fluctuate in sign, which
                ! means that they would threaten to spoil the calculation if they were included.
                ! In theory, they should also have zero frequencies, as they represent translations
                ! and rotations of the cluster.  We should thus achieve a good approximation to the
                ! ZPE if the lowest six eigenvalues are omitted from the calculation.

                WRITE(99,*) 'frequencies in 1/cm = '
                WRITE(99,*) sqrt(abs(evalues)) * aucm
                WRITE(99,*) 'Vibrational zero-point energy =', sum(sqrt(evalues(6:)))/2, ' hartree'
                WRITE(99,*) '                              =', (sum(sqrt(evalues(6:)))/2)*aucm, ' 1/cm'
                WRITE(*,*) 'frequencies in 1/cm = '
                WRITE(*,*) sqrt(abs(evalues)) * aucm
                WRITE(*,*) 'Vibrational zero-point energy =', sum(sqrt(evalues(6:)))/2, ' hartree'
                WRITE(*,*) '                              =', (sum(sqrt(evalues(6:)))/2)*aucm, ' 1/cm'
        end subroutine

!=============================================================================================================================

        subroutine do_hessian_multiple(evalues)
                INTEGER :: ldz = 2                                      ! No. of evectors to include in evectors matrix
                double precision,DIMENSION(:,:),ALLOCATABLE :: hessian
                double precision,DIMENSION(:,:),ALLOCATABLE :: evectors	! Computed eigenvectors of matrix
                double precision,DIMENSION(:),ALLOCATABLE,intent(out) :: evalues	! Computed eigenvalues of matrix

                ! =======================================================================
                ! CALCULATE HESSIAN
                ! =======================================================================

                IF (myrank == root) THEN                        ! To avoid duplication in non-root processes
                        WRITE(99,'(A)') '-----------------------------------------------------------------------'
                        WRITE(99,'(A)') 'Calculating Hessian matrix of second derivatives'
                        WRITE(99,'(A)') '-----------------------------------------------------------------------'
                        WRITE(*,'(A)') '-----------------------------------------------------------------------'
                        WRITE(*,'(A)') 'Calculating Hessian matrix of second derivatives'
                        WRITE(*,'(A)') '-----------------------------------------------------------------------'
                        WRITE(*,'(A,ES10.3,A)') 'HESSEPS for Hessian calculation set to ', hesseps, ' bohr'
                        WRITE(99,'(A,ES10.3,A)') 'HESSEPS for Hessian calculation set to ', hesseps, ' bohr'

                        ALLOCATE( hessian(((d*p)+1),d*p*n) )        ! Create the blank array for the Hessian matrix of second derivatives
                        hesst = MPI_Wtime()
                        CALL hesscalc ( x, hessian, .true. ) ! setting the last parameter to .false. restores the 
                                                                ! old behaviour of the serial routine
                        hesst = MPI_Wtime() - hesst
                else
                        call hesscalc_child( x )     !special routine for child processes 
                                                        ! (doesn't need an allocated hessian matrix)
                        return
                endif
                !TODO: Diagonalisation is not parallelised!!
                ! This is only done by root process

                ! =======================================================================
                ! DIAGONALIZE HESSIAN
                ! =======================================================================

                WRITE(99,'(A)') '-----------------------------------------------------------------------'
                WRITE(99,'(A)') 'Diagonalizing Hessian matrix of second derivatives'
                WRITE(99,'(A)') '-----------------------------------------------------------------------'
                WRITE(*,'(A)') 'Diagonalizing Hessian matrix of second derivatives'
                WRITE(*,'(A)') '-----------------------------------------------------------------------'

                hessDiagT = MPI_Wtime()
                ALLOCATE( evectors(ldz,d*p*n) )		! Create the blank array for the eigenvectors of the Hessian
                ALLOCATE( evalues(d*p*n) )		! Create the blank array for the eigenvalues of the Hessian
                !CALL hessout ( d, p, n, hessian, hessianoutput )	! To write the Hessian to a file, for debugging purposes
                CALL diagband ( hessian, evalues, evectors, ldz )

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !CALL externalevals	! To override calculated evals for testing purposes
                !evalues = extevals
                !WRITE(99,'(A)') '-----------------------------------------------------------------------'
                !WRITE(99,'(A)') 'WARNING - EIGENVALUES OVERRIDDEN - INPUT IS FROM FILE CALLED EIGVALS'
                !WRITE(99,'(A)') '-----------------------------------------------------------------------'
                !WRITE(*,'(A)') 'WARNING - EIGENVALUES OVERRIDDEN - INPUT IS FROM FILE CALLED EIGVALS'
                !WRITE(*,'(A)') '-----------------------------------------------------------------------'
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                hessDiagT = MPI_Wtime() - hessDiagT

                WRITE(99,'(A)')
                WRITE(99,*) 'evalues = '
                WRITE(99,*) evalues
                WRITE(99,'(A)')
                !WRITE(*,'(A)')
                !WRITE(*,*) 'evalues(:12) = '
                !WRITE(*,*) evalues(:12)
                !WRITE(*,'(A)')
        end subroutine

!=============================================================================================================================

        subroutine do_splitting_multiple(evalues,skink)
                double precision,DIMENSION(:),ALLOCATABLE,intent(in) :: evalues    ! Computed eigenvalues of matrix
                double precision, allocatable,dimension(:) :: zeroedevalues
                double precision, allocatable,dimension(:) :: evalsforphi
                double precision,intent(in) :: skink       ! Computed eigenvalues of matrix
                double precision :: logdet     ! To contain the log of the determinant (sum of logs of non-zero evalues)

                double precision,dimension(:),allocatable :: evalues_collapsed !contains eigenvalues of collapsed linear polymer
                double precision,dimension(:),allocatable :: buffer !temporary buffer for above data
                double precision :: logdet_collapsed
                double precision :: sinterm
                double precision :: phi, h
                integer :: i,j,offs

                !No real need to paralellise this:
                if ( myrank .ne. root) return

                WRITE(99,'(A)') ''
                WRITE(99,'(A)') '######################################################################'
                WRITE(99,'(A)') "#---                 SPLITTING CALCULATION                        ---#"
                WRITE(99,'(A)') '######################################################################'
                WRITE(*,'(A)') ''
                WRITE(*,'(A)') '######################################################################'
                WRITE(*,'(A)') "#---                 SPLITTING CALCULATION                        ---#"
                WRITE(*,'(A)') '######################################################################'

                ! =======================================================================
                ! CALCULATE LOGDET
                ! =======================================================================
                ! "...the pseudo-determinant is the product of all non-zero eigenvalues of a square matrix.
                ! It coincides with the regular determinant when the matrix is non-singular." (Wikipedia)
                !
                ! For our case (1 kink) there is only 1 nearly-zero eigenvalue ( evalues(1) )
                ! We do *not* remove the 6 trans/rot eigenvalues ( evalues(2:7) )
                ! => interested in:
! DEBUG
	        logdet = SUM( dLOG(evalues(2:)) )
!                logdet = SUM( dLOG(evalues(8:)) )

!                logdet = -38363.2445548403
!                ALLOCATE( evalsforphi((d*p*n)-1) )		! Create the blank array for the eigenvalues of the Hessian
!                evalsforphi(1:6) = evalues(1:6)
!                evalsforphi(7:(d*p*n)-1) = evalues(8:(d*p*n))
!                logdet = SUM( dLOG(evalsforphi) )
!                WRITE(*,*) 'evalues(1:10)', evalues(1:10)
!                WRITE(*,*) 'evalsforphi(1:10)', evalsforphi(1:10)
! DEBUG
!                WRITE(*,*) 'OVERRIDING USUAL PHI CALCULATION!'
!                WRITE(99,*) 'OVERRIDING USUAL PHI CALCULATION!'

                ! UNITS OF LOGDET:
                ! logdet = prod(Jevalues(2:)), i.e. one eigenvalue has been removed from the product
                ! The units of the eigenvalues are found as follows:
                !  - The Planck relation tells us that E = hbar * omega
                !  - This rearranges to omega = E / hbar
                !  - The eigenvalues are equal to omega squared, hence they have units (hartree/hbar)^2
                !  - There are (d*p*n)-1 eigenvalues in logdet
                !  - Thus, the units are [(hartree/hbar)^2]^[(d*p*n)-1] = (hartree/hbar)^2[(d*p*n)-1]

                WRITE(*,*) 'Pseudo-determinant of instanton found via SUM( LOG( evalues(2:) ) )'
                WRITE(*,'(A,E15.8,A,I5)') '       log(det) = ', logdet, ' (hartree/hbar) ^', 2*((cdata%dof*n)-1)
                WRITE(99,*) 'Pseudo-determinant of instanton found via SUM( LOG( evalues(2:) ) )'
                WRITE(99,'(A,E15.8,A,I5)') '       log(det) = ', logdet, ' (hartree/hbar) ^', 2*((cdata%dof*n)-1)

                ! =======================================================================
                ! CALCULATE PHI
                ! =======================================================================
                allocate(evalues_collapsed(cdata%dof*n))
                allocate(buffer(cdata%dof*n))
                evalues_collapsed =0d0
                buffer=0d0

                allocate(zeroedevalues(cdata%dof))
                zeroedevalues = cdata%evalues
! DEBUG
                zeroedevalues(:6) = 0d0	! Set lowest 6 eigenvalues to zero
                ! This corrects for errors in the calculation and diagonalization of the
                ! single-bead Hessian, as in principle the six lowest eigenvalues should
                ! be exactly zero (translations and rotations), whereas in practice they
                ! end up as very small (but non-zero) values, including negative numbers.

                ! For this routine see eq. (5.22) of Jeremy's thesis
                ! Essentially we calculate all sine terms and take all zeroedevalues
                ! and form all possible combinations and sort them
                do i=1,n        !loop over beads
                        sinterm=( 2d0/(betan*hbar) * dsin(pi*i/(2*n+2)) )**2
                        offs=(i-1)*cdata%dof
                        buffer(offs+1:offs+cdata%dof)=zeroedevalues+sinterm
                        ! Over the range we use the sin it is monotonic
                        ! since the evalues were saved (and read in) in ascending order
                        ! we can hence assume that the section of the buffer
                        ! just filled is sorted ascendingly.
                        ! 
                        ! By construction evalues_collapsed is already sorted,
                        ! so we can simply merge these two:
                        call dmerge(offs,evalues_collapsed,cdata%dof,buffer(offs+1),offs+cdata%dof,buffer)
                        ! Note: buffer(offs) sets the array up s.t. the first
                        ! element taken by dmerge for the second array is the
                        ! one at offs+1
                        evalues_collapsed = buffer
                enddo !i
                ! TODO: This can be even more sped up if we don't use a linear
                ! merge (chunk merged into sorted array one after another) but
                ! if we used a logarithmic merge (chunks from pairs, two smaller
                ! pairs are merged; the larger, merged chunks form pairs again
                ! and are merged pairwise, ... until the full array has been
                ! merged.
                deallocate(buffer)
               
                if (debug) then
                        write(99,*)  
                        write(99,*) "Evalues of collapsed linpol:"
                        do i=1,cdata%dof*n
                                write(99,*) evalues_collapsed(i)
                        enddo
                end if

                logdet_collapsed = SUM( dLOG(evalues_collapsed(:)) )
! DEBUG
!                logdet_collapsed = SUM( dLOG(evalues_collapsed(7:)) )

                if (debug) then
                        write(99,*)
                        write(99,*) "logdet linpol collapsed:"
                        write(99,*) logdet_collapsed
                end if

                ! UNITS OF LOGDET_COLLAPSED:
                ! logdet_collapsed = prod(J0evalues(:)), i.e. all eigenvalues are present in the product
                ! The units of the eigenvalues are found as follows:
                !  - The Planck relation tells us that E = hbar * omega
                !  - This rearranges to omega = E / hbar
                !  - The eigenvalues are equal to omega squared, hence they have units (hartree/hbar)^2
                !  - There are d*p*n eigenvalues in logdet_collapsed
                !  - Thus, the units are [(hartree/hbar)^2]^(d*p*n) = (hartree/hbar)^2(d*p*n)

                WRITE(*,*) 'Pseudo-determinant of collapsed LP found via SUM( LOG( evalues(:) ) )'
                WRITE(*,'(A,E15.8,A,I5)') '       log(det) of collapsed polymer = ', &
                        logdet_collapsed, ' (hartree/hbar)^', 2*cdata%dof*n
                WRITE(99,*) 'Pseudo-determinant of collapsed LP found via SUM( LOG( evalues(:) ) )'
                WRITE(99,'(A,E15.8,A,I5)') '       log(det) of collapsed polymer = ', &
                        logdet_collapsed, ' (hartree/hbar)^', 2*cdata%dof*n

                ! UNITS OF PHI:
                ! phi = sqrt((det J')/(det J0)) = sqrt( prod(Jevalues(2:)) / prod(J0evalues(:)) )
                ! Jevalues(2:) contains one fewer eigenvalues than J0evalues(:).  Hence, we divide their units to find:
                ! (hartree/hbar)^2[(d*p*n)-1] / (hartree/hbar)^2(d*p*n) = (hartree/hbar)^-2 = (hbar/hartree)^2
                ! We complete the derivation by taking the root, to obtain units of hbar/hartree

                phi = dexp ( (logdet - logdet_collapsed) / 2.0 )
                write(*,*) 'Calculating ratio of tunnelling and non-tunnelling determinants ...'
                write(99,*) 'Calculating ratio of tunnelling and non-tunnelling determinants ...'
                WRITE(*,*) '       phi = ', phi, ' hbar/hartree'
                WRITE(99,*) '       phi = ', phi, ' hbar/hartree'

                ! =======================================================================
                ! CALCULATE h
                ! =======================================================================
                ! Here we calculate an estimate of h for the values of Skink, phi
                ! assuming we are in the limiting case of (beta -> Infinity)
                ! 
                ! From Dimer/Trimer Paper eqns (17, 8) we can derive the expression
                !      h = -hbar/Phi * sqrt(Skink/(2*pi*hbar)) * exp(-Skink/hbar)
                ! where hbar = 1 (Atomic units)

                h = -hbar/phi * dsqrt(skink/(2d0*pi*hbar)) * dexp(-skink/hbar)

                ! units of h = hbar/(hbar/hartree) * sqrt(hbar/hbar) * exp(hbar/hbar)
                !            = hartree

                write(*,*) 'Calculating kink path weight ...'
                write(99,*) 'Calculating kink path weight ...'
                WRITE(*,*) '      h = ', h, ' hartree'
                WRITE(*,*) '        = ', h*aucm, ' 1/cm'
                WRITE(99,*) '      h = ', h, ' hartree'
                WRITE(99,*) '        = ', h*aucm, ' 1/cm'
        end subroutine

!=============================================================================================================================

        subroutine do_end(successfulCompletion) ! display runtime and end program:
                logical :: successfulCompletion ! Has code finished successfully
                character(len=46) :: dateString

                ! =======================================================================
                ! STUFF DONE BY NON-ROOT
                ! =======================================================================

                if ( myrank .ne. root) then
                        !just finalise mpi and end
                        CALL MPI_Finalize(error)
                        stop
                end if
                !only done by root

                ! =======================================================================
                ! CLOSING TEXT
                ! =======================================================================
                IF (successfulCompletion) THEN
                        WRITE(99,'(A)')
                        WRITE(99,'(A)') '======================================================================='
                        WRITE(99,'(A)') 'Code complete'
                        WRITE(99,'(A)') '======================================================================='
                        WRITE(99,'(A)')
                        WRITE(*,'(A)') '======================================================================='
                        WRITE(*,'(A)') 'Code complete'
                        WRITE(*,'(A)') '======================================================================='
                        WRITE(*,'(A)')
                else
                        WRITE(99,'(A)')
                        WRITE(99,'(A)') '======================================================================='
                        WRITE(99,'(A)') 'Code exited without completing successfully'
                        WRITE(99,'(A)') '======================================================================='
                        WRITE(*,'(A)')
                        WRITE(*,'(A)') '======================================================================='
                        WRITE(*,'(A)') 'Code exited without completing successfully'
                        WRITE(*,'(A)') '======================================================================='
                        WRITE(*,'(A)')

                end if

                ! =======================================================================
                ! RECORD MPI END TIME, DISPLAY TOTAL MPI TIME, FINALIZE MPI
                ! =======================================================================

                startt = (MPI_Wtime()-startt)

                write(* ,'(A)') "Runtime (per process) was:"
                write(99,'(A)') "Runtime (per process) was:"
                dateString = dateToString(mint)
                write(* ,'(3x,2A)') "Minimisation (LBFGS):    ", dateString
                write(99,'(3x,2A)') "Minimisation (LBFGS):    ", dateString
                dateString = dateToString(hesst)
                write(* ,'(3x,2A)') "Hessian calculation:     ", dateString
                write(99,'(3x,2A)') "Hessian calculation:     ", dateString
                dateString = dateToString(hessDiagT)
                write(* ,'(3x,2A)') "Hessian diagonalisation: ", dateString
                write(99,'(3x,2A)') "Hessian diagonalisation: ", dateString
                write(* ,'(3x,A)')  "---------------------------------------------------------------------"
                write(99,'(3x,A)')  "---------------------------------------------------------------------"
                dateString = dateToString(startt)
                write(* ,'(3x,2A)') "Overall runtime:         ", dateString
                write(99,'(3x,2A)') "Overall runtime:         ", dateString
                WRITE(*,*)
                WRITE(99,*)

                CALL MPI_Finalize(error)
                stop
        end subroutine
END PROGRAM linpolexe
