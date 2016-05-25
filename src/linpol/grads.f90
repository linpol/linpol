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
MODULE grads
USE parameters
USE funcs
USE prep
use pes_interface
implicit none
!
! =======================================================================
! PURPOSE OF THIS MODULE
! =======================================================================
! To contain the subroutines to calculate the gradients of the internal
! and external ring-polymer potentials, plus the gradient of the overall
! total potential.
! 
! Contents:
!     intgradscaled
!     extgrad
!     totgradx
!     onebeadgradx
!
! Record of revisions:
!
!  Date           Person         Change implemented
!  ----           ------         ------------------
!  21 Nov 2011    A. Reid        Development of module
!  29 Nov 2011    A. Reid        Conversion to use variables from modules
!  1 June 2012    A. Reid        Redundant subroutines deleted
!
! =======================================================================

CONTAINS

! ===========================================================================================================

    SUBROUTINE intgradscaled ( xinput, internalgrad )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To calculate the 'internal' part of the ring-polymer potential gradient
    ! 
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  9 Jan 2012     A. Reid        Development of subroutine
    !  9 May 2012     A. Reid        Fixed ends option commented out
    !
    ! =======================================================================

    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================

    ! Calling variables:
    double precision,DIMENSION(d,p,n),INTENT(IN) :: xinput            ! Linear polymer geometry
    double precision,DIMENSION(d*p*n),INTENT(OUT) :: internalgrad     ! Internal linear-polymer gradient

    ! Local variables:
    double precision,DIMENSION(d,p,n) :: x                            ! Linear polymer geometry
    INTEGER :: dim                                                 ! Index of DO loops
    INTEGER :: atom                                                ! Index of sub DO loops
    INTEGER :: i                                                   ! Index of beads
    double precision :: prefactor                                     ! 
    double precision :: grad                                          ! 

    x = xinput

    prefactor = 1.0D0 / ( (betan*hbar)**2 )                    ! Pre-factor for the first derivative - note lack of 2 in denominator

    beadloop: DO i=1,n 
        atomloop: DO atom=1,p
            dimloop: DO dim=1,d

                !When i = 1:
                IF (i==1) THEN
                    SELECT CASE(endchoice)
                        CASE(1)  ! fixed ends
                            grad = prefactor * (masses(atom)) * (  2.0d0*x(dim,atom,1)  -  x(dim,atom,2)  -  xA(dim,atom))
!                        CASE(2)  ! Unfixed ends
!                            !prefactor*(x_l,1 - x_l,2)
!                            grad = prefactor * (masses(atom)) * (  x(dim,atom,1)  -  x(dim,atom,2)  )
                        CASE DEFAULT  ! endchoice not specified - defaults to unfixed ends
                            !prefactor*(x_l,1 - x_l,2)
                            grad = prefactor * (masses(atom)) * (  x(dim,atom,1)  -  x(dim,atom,2)  )
                    END SELECT
                ENDIF

                !When i = L:
                IF (i==n) THEN

                    SELECT CASE(endchoice)
                        CASE(1)  ! fixed ends
                            grad = prefactor * (masses(atom)) * (  2.0d0*x(dim,atom,n)  -  x(dim,atom,n-1)  -  xB(dim,atom))
!                        CASE(2)  ! Unfixed ends
!                            !NEW:prefactor*(x_l,L - x_l,L-1)
!                            grad = prefactor * (masses(atom)) * (  x(dim,atom,n)  -  x(dim,atom,n-1)  )
                        CASE DEFAULT  ! endchoice not specified - defaults to unfixed ends
                            !NEW:prefactor*(x_l,L - x_l,L-1)
                            grad = prefactor * (masses(atom)) * (  x(dim,atom,n)  -  x(dim,atom,n-1)  )
                    END SELECT

                ENDIF
                IF ((i.LT.n).AND.(i.GT.1)) THEN  ! The same regardless of whether the ends are
                    !When i =/ 1,L:
                    !prefactor*(2x_l,i - x_l,i-1 - x_l,i+1)
                    grad = prefactor * (masses(atom)) * (  ( 2.0d0*x(dim,atom,i) )  -  x(dim,atom,i-1) -   x(dim,atom,i+1)  )
                ENDIF
                internalgrad( calcpos ( d, p, dim, atom, i ) ) = grad
            END DO dimloop
        END DO atomloop
    END DO beadloop

    END SUBROUTINE intgradscaled

! ===========================================================================================================

    SUBROUTINE extgrad ( xin, externalgrad )
    USE mpi
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To calculate the gradient of the 'external' part of the linear-polymer
    ! potential.
    ! 
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  1 Nov 2011     A. Reid        Development of subroutine
    !  24 May 2012    A. Reid        Parallelization with MPI
    !
    ! =======================================================================

    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================

    ! Calling variables:
    double precision,DIMENSION(d,p,n),INTENT(IN) :: xin		! Ring polymer geometry
    double precision,DIMENSION(d*p*n),INTENT(OUT) :: externalgrad	! External linear-polymer gradient

    ! Local variables:
    double precision,DIMENSION(:),allocatable :: beadgrad		! Gradient array for single bead
    INTEGER :: i						! Index of time slices
    INTEGER :: ind						! Index of DO loops

    ! MPI variables:
    double precision,DIMENSION(:,:),allocatable :: beadsendbuf	! Linear polymer geometry split into n beads, to be scattered to processes
    double precision,DIMENSION(:),allocatable :: beadrecvbuf	! Bead geometries, to be received by one of nproc processes
    double precision,DIMENSION(:),allocatable :: beadgeom		! Linear polymer geometry for calculation of potential
    double precision,DIMENSION(:),allocatable :: gradsendbuf  	! Gradient of all beads in one of nprocs processes
    double precision,DIMENSION(:),allocatable :: gradrecvbuf	! Gradient for all linear polymer beads, to be gathered from nprocs processes
    double precision :: potential				! Sum of single bead potentials, non-zero except in the root process
    INTEGER :: error, nprocs, namelen, myrank

    ! =======================================================================
    ! ALLOCATE SOME ARRAYS
    ! =======================================================================
    allocate(beadgeom(d*p))
    allocate(beadsendbuf(d*p,n))
    allocate(gradrecvbuf(d*p*n))
    allocate(beadgrad(d*p))

    ! =======================================================================
    ! GENERATE MPI DATA
    ! =======================================================================
    CALL MPI_Comm_size ( MPI_COMM_WORLD, nprocs, error)
    CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)

    ! =======================================================================
    ! PROCEED WITH MESSAGE-PASSING
    ! =======================================================================

    ! Initial setup of the linear polymer bead array within the root process only
    initialsetupmulti: IF (myrank == root) THEN
        CALL threetoone ( d, p, n, xin, externalgrad )
        DO i = 1, n
            beadsendbuf(1:(d*p),i) = externalgrad( ((i-1)*(d*p))+1 : (i*(d*p)) )
        END DO
    END IF initialsetupmulti


    ALLOCATE(beadrecvbuf(d*p*n/nprocs))
    ALLOCATE(gradsendbuf(d*p*n/nprocs))

    ! Initial setup of the x array is complete, so we now scatter the beads to the different processes
    CALL MPI_Scatter(beadsendbuf, d*p*n/nprocs, MPI_DOUBLE_PRECISION, &
         beadrecvbuf, d*p*n/nprocs, MPI_DOUBLE_PRECISION, root, &
         MPI_COMM_WORLD, error)
           
    DO i = 1, n/nprocs
        beadgeom = beadrecvbuf( ((i-1)*(d*p))+1 : (i*(d*p)) )
        call pes_obj%grad(beadgeom,gradeps,beadgrad)
        !beadgrad=grad(beadgeom,gradeps)			! See pes_init for function
        gradsendbuf( ((i-1)*p*d)+1 : i*p*d ) = beadgrad	! Running total of gradient for all beads in this process
    END DO

    ! Calculation of bead gradients is complete, so we now gather them back from the different processes
    CALL MPI_Gather(gradsendbuf, d*p*n/nprocs, MPI_DOUBLE_PRECISION, &
         gradrecvbuf, d*p*n/nprocs, MPI_DOUBLE_PRECISION, root, &
         MPI_COMM_WORLD, error)

    ! The root process now has all the bead gradients, so we sum them to find the overall gradient
    IF (myrank == root) THEN
        externalgrad = gradrecvbuf
    ENDIF

    ! We now broadcast externalgrad to all processes
    CALL MPI_Bcast(externalgrad, d*p*n, MPI_DOUBLE_PRECISION, &
         root, MPI_COMM_WORLD, error)

    deallocate(gradsendbuf,beadrecvbuf, beadgeom,gradrecvbuf,beadgrad,beadsendbuf)
    END SUBROUTINE extgrad

! ===========================================================================================================

    SUBROUTINE totgradx ( xin, totalgrad )
    !TODO: There is a bug in here giving a segmentation fault for large systems

    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To calculate the total gradient of the ring-polymer potential using the
    ! subroutines for the internal and external parts.  Takes input in format
    ! compatible with Marko's L-BFGS code.
    ! 
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  23 Nov 2011    A. Reid        Development of subroutine
    !
    ! =======================================================================

    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================

    ! Calling variables:
    double precision,DIMENSION(d*p*n),INTENT(IN) :: xin                 ! Ring polymer geometry
    double precision,DIMENSION(d*p*n),INTENT(OUT) :: totalgrad      ! Total linear-polymer gradient

    ! Local variables:
    double precision,DIMENSION(d,p,n) :: x              ! Linear polymer geometry    
    double precision,DIMENSION(d,p,n) :: xscaled        ! Scaled linear polymer geometry
    double precision,DIMENSION(d*p*n) :: internalgrad   ! Internal linear-polymer gradient
    double precision,DIMENSION(d*p*n) :: externalgrad   ! External linear-polymer gradient

    ! =======================================================================
    ! CONVERTING FROM 1D ARRAY TO 3D ARRAY
    ! =======================================================================

    CALL onetothree ( d, p, n, xin, x )

    ! =======================================================================
    ! GRADIENT OF LP POTENTIAL ('INTERNAL' PART)
    ! =======================================================================

    CALL intgradscaled ( x, internalgrad )

    ! =======================================================================
    ! GRADIENT OF LP POTENTIAL ('EXTERNAL' PART)
    ! =======================================================================

    CALL extgrad ( x, externalgrad )

    ! =======================================================================
    ! GRADIENT OF LP POTENTIAL (TOTAL)
    ! =======================================================================

    totalgrad = internalgrad + externalgrad

    END SUBROUTINE totgradx

! ===========================================================================================================

!     SUBROUTINE onebeadgradx ( xonebead, onebeadgradient )
!     USE pes_shell
!     !
!     ! =======================================================================
!     ! PURPOSE OF THIS SUBROUTINE
!     ! =======================================================================
!     ! To calculate the gradient of the 'external' potential for a single
!     ! bead.
!     ! 
!     !
!     ! Record of revisions:
!     !
!     !  Date           Person         Change implemented
!     !  ----           ------         ------------------
!     !  5 Dec 2011     A. Reid        Development of subroutine
!     !
!     ! =======================================================================
! 
!     ! =======================================================================
!     ! DECLARATION OF VARIABLES
!     ! =======================================================================
! 
!     ! Calling variables:
!     REAL(kind=dp),DIMENSION(d*p),INTENT(IN) :: xonebead                 ! Ring polymer geometry
!     REAL(kind=dp),DIMENSION(d*p),INTENT(OUT) :: onebeadgradient     ! External linear-polymer gradient
! 
!     ! Local variables:
!     REAL(kind=dp),DIMENSION(d*p) :: xcalc 
! 
!     xcalc=xonebead
!     onebeadgradient=grad(xcalc,gradeps)                                 ! See pes_shell for function
! 
!     END SUBROUTINE onebeadgradx

! ===========================================================================================================

    function onebeadgnorm(x) result(norm)
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To calculate the gradient modulus of the 'external' potential for a single
    ! bead.
    ! 
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  20 Jan 2013    M. Herbst      Development of subroutine
    !
    ! =======================================================================

    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================

    ! Calling variables:
    double precision,DIMENSION(:),INTENT(IN) :: x                ! Ring polymer geometry

    ! Local variables:
    double precision,DIMENSION(size(x)) :: gradient      ! External linear-polymer gradient
    double precision :: norm
    integer :: i

    norm=0.0
    call pes_obj%grad(x,gradeps,gradient)

      !WRITE(*,*) 'DEBUG - Gradient vector:'
      !WRITE(*,*) gradient

    do i=1,size(x)
        norm = norm + (gradient(i))**2
    end do
    norm=sqrt(norm)
    END function

END MODULE grads
