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
MODULE min
USE parameters						! Module containing parameters to be pulled into all other modules
USE prep						! Linear polymer preparation module
USE pots						! Potential calculation module
USE grads						! Gradient calculation module
use opt_interface                                       ! interface for optimisations
use opt_selected
!
! =======================================================================
! PURPOSE OF THIS MODULE
! =======================================================================
! To contain the functions and subroutines to accomplish miscellaneous
! tasks.  Had to be set up as a separate module because it relies upon
! subroutines and data from the other modules.
! 
! Contents:
!     minsingle
!     minlinpol
!
! Record of revisions:
!
!  Date           Person         Change implemented
!  ----           ------         ------------------
!  22 Dec 2011    A. Reid        Development of module
!  9 May 2012     A. Reid        Removed old commented-out code
!  9 May 2012     A. Reid        Format descriptors tidied up and IF comments removed
!  9 May 2012     A. Reid        'minstartend' subroutine deleted
!
! =======================================================================

IMPLICIT NONE
SAVE                                       ! Ensures all declared data values will be preserved

! =======================================================================
! DECLARATION OF GLOBAL VARIABLES AVAILABLE TO OTHER PROGRAM UNITS
! =======================================================================

    double precision,DIMENSION(:),ALLOCATABLE :: xsinglemin	! Minimum energy starting geometry
    double precision :: totalpot					! Total linear-polymer potential
    double precision :: singlepot					! Single bead potential

CONTAINS

    SUBROUTINE minsingle ( outputfile )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To find the minimum energy geometries of any configurations of a water
    ! cluster.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  18 Mar 2012    A. Reid        Development of subroutine
    !  9 May 2012     A. Reid        Name changed from 'minstart' to 'minsingle'
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Calling variables:
    CHARACTER(len=64),INTENT(IN) :: outputfile
    
    ! Local variables:
    double precision,DIMENSION(:),ALLOCATABLE :: singlegrad	! Overall starting geometry gradient


    ! ALLOCATING ARRAYS

    ALLOCATE( singlegrad( d*p ) )			! Create the blank array for the overall starting geometry gradient
    ALLOCATE( xsinglemin( d*p ) )			! Create the blank minimum energy starting geometry

    ! FINDING MINIMUM ENERGY STARTING GEOMETRY
    WRITE(*,'(A,ES10.3,A)') 'EPS for optimization set to ', starteps, ' hartree/bohr'
    WRITE(*,'(A,ES10.3,A)') 'GRADEPS for finite-difference gradient set to ', gradeps, ' bohr'
    WRITE(99,*)
    WRITE(99,'(A,ES10.3,A)') 'EPS for optimization set to ', starteps, ' hartree/bohr'
    WRITE(99,'(A,ES10.3,A)') 'GRADEPS for finite-difference gradient set to ', gradeps, ' bohr'
    WRITE(*,*)

    !##########################################
    !#-- transform beads to pes coordinates --#
    !##########################################
    call pes_obj%atomvec_to_pes(masses)
    call pes_obj%vec_to_pes(xsingle)
    !call pes_obj%vec_to_pes(atomtags)
    transformed_to_pot = .true.
    !##########################################

    xsinglemin = RESHAPE(xsingle, (/ d*p /))
    singleBeadMode = .true. !TODO: Temporary
    !msave: number of steps to save in the lbfgs iteration
    call opt_obj%optimise( xsinglemin, singlepot, singlegrad, msave=3, eps=starteps, maxiterations=lbfgs_maxiter )
    xsingle = RESHAPE(xsinglemin, (/d,p/))

    !############################################
    !#-- transform beads from pes coordinates --#
    !############################################
    call pes_obj%vec_from_pes(xsingle)
    call pes_obj%atomvec_from_pes(masses)
    !call pes_obj%vec_from_pes(atomtags)
    transformed_to_pot = .false.
    !############################################

    CALL startendxyzout ( xsingle, d, p, beta, atomtags, outputfile )
    END SUBROUTINE minsingle

! ===========================================================================================================

    SUBROUTINE minlinpol ( outputfile )
    USE mpi
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To find the minimum energy geometries of any configurations of a water
    ! cluster.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  21 May 2012    A. Reid        Development of subroutine from 'minsingle'
    !  1 June 2012    A. Reid        totpot/totgrad replaced by totpotx/totgradx
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Calling variables:
    CHARACTER(len=64),INTENT(IN) :: outputfile
    
    ! Local variables:

    double precision,DIMENSION(:),ALLOCATABLE :: totalgrad	! Total linear-polymer gradient
    double precision,DIMENSION(:),ALLOCATABLE :: xminone	! Unscaled minimum-energy linear-polymer geometry (1D array)
    integer :: i

    INTEGER :: error, nprocs, namelen, myrank

    ! =======================================================================
    ! GENERATE MPI DATA
    ! =======================================================================
    CALL MPI_Comm_size ( MPI_COMM_WORLD, nprocs, error)
    CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)


    ! ALLOCATING ARRAYS

    ALLOCATE( totalgrad( d*p*n ) )                              ! Create the blank array for the total linear-polymer gradient
    ALLOCATE( xminone( d*p*n ) )                                ! Create the blank 1-D array to hold the minimum-energy geometry

    IF (myrank == root) THEN				! To avoid duplication in non-root processes
        WRITE(*,*)
        WRITE(*,'(A,ES10.3,A)') 'EPS for optimization set to ', choseneps, ' hartree/bohr'
        WRITE(*,'(A,ES10.3,A)') 'GRADEPS for finite-difference gradient set to ', gradeps, ' bohr'
        WRITE(99,*)
        WRITE(99,'(A,ES10.3,A)') 'EPS for optimization set to ', choseneps, ' hartree/bohr'
        WRITE(99,'(A,ES10.3,A)') 'GRADEPS for finite-difference gradient set to ', gradeps, ' bohr'
    END IF

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

    xminone = RESHAPE(x, (/ d*p*n /))
    CALL opt_obj%optimise ( xminone, totalpot, totalgrad, msave=3, eps=choseneps, maxiterations=lbfgs_maxiter )
    x = RESHAPE(xminone, (/ d, p, n /))

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

    CALL xyzout ( x , n , d , p, beta, atomtags, outputfile )

    deallocate( xminone, totalgrad)
    END SUBROUTINE minlinpol

! ===========================================================================================================

END MODULE min
