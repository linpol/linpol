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
MODULE pots
USE parameters
USE funcs
USE prep
use pes_interface
USE clusterdatamod
implicit none
!
! =======================================================================
! PURPOSE OF THIS MODULE
! =======================================================================
! To contain the subroutines to calculate the internal and external
! ring-polymer potentials, plus the overall total potential.
! 
! Contents:
!     intpotscaled
!     extpot
!     totpotx
!     onebeadpotx
!     pathway
!
! Record of revisions:
!
!  Date           Person         Change implemented
!  ----           ------         ------------------
!  21 Nov 2011    A. Reid        Development of module
!  29 Nov 2011    A. Reid        IF statements added to govern comments
!  29 Nov 2011    A. Reid        Conversion to use variables from modules
!  1 June 2012    A. Reid        Redundant subroutines deleted
!
! =======================================================================

CONTAINS

! ===========================================================================================================

    SUBROUTINE intpotscaled ( xinput, internal )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To calculate the 'internal' part of the ring-polymer potential.  Based
    ! on the old 'intpot' subroutine, but takes an unscaled array rather than
    ! a scaled one.
    ! 
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  6 Jan 2012     A. Reid        Development of subroutine
    !  9 May 2012     A. Reid        Fixed ends option commented out
    !
    ! =======================================================================

    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================

    ! Calling variables:
    double precision,DIMENSION(d,p,n),INTENT(IN) :: xinput      ! Linear polymer geometry
    double precision,INTENT(OUT) :: internal                    ! Internal linear polymer potential

    ! Local variables:
    double precision,DIMENSION(d,p,n) :: x                      ! Linear polymer geometry
    INTEGER :: dim                                           ! Index of DO loops
    INTEGER :: atom                                          ! Index of sub DO loops
    INTEGER :: bead                                          ! Index of sub-sub DO loops
    double precision :: massroot                                ! Square root of the mass of the degree of freedom being dealt with
    double precision :: fixedterms                              ! Fixed terms from 'internal' sum
    double precision :: variableterms                           ! Variable terms from 'internal' sum
    double precision :: calcf                                   ! Temporary value to find 'fixedterms'
    double precision :: calcv                                   ! Temporary value to find 'variableterms'
    double precision :: prefactor                               ! 1/(2*(betan*hbar)**2)

    ! Setting terms to zero before start of calculation run
    fixedterms = 0.0d0
    variableterms = 0.0d0
    calcf = 0.0d0
    calcv = 0.0d0

    prefactor = 1.0D0 / ( 2.0D0 * ((betan*hbar)**2) )                    ! Pre-factor for the linear-polymer potential

    x = xinput

    ! 'FIXED' TERMS (CAN ALSO BE UNFIXED, I.E. TERM SET TO ZERO)

    SELECT CASE(endchoice)  ! 'endchoice' variable is set in the parameters.f90 file
        CASE(1) !fixed ends
            fixedterms = sum(masses(:) * sum((x(:,:,1) - xA(:,:))**2 + (x(:,:,n) - xB(:,:))**2,1))
        CASE DEFAULT  ! endchoice not specified - defaults to unfixed ends
            fixedterms = 0.0d0
    END SELECT


    ! VARIABLE TERMS

    overdimsv: DO dim = 1, d
        overatomsv: DO atom = 1, p

            variablecalc: DO bead = 1, n-1, 1                        ! DO loop from 1 to n-1 in increments of 1
                calcv = (masses(atom))   *   (  ( x(dim,atom,bead+1)-x(dim,atom,bead) )**2  )
                variableterms = variableterms + calcv
            END DO variablecalc

        END DO overatomsv
    END DO overdimsv

    internal = prefactor * (fixedterms + variableterms)

    END SUBROUTINE intpotscaled

! ===========================================================================================================

    SUBROUTINE extpot ( xin, potsum )
    USE mpi
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To calculate the 'external' part of the ring-polymer potential.  Note
    ! that the input geometry must have already been converted from Angstroms
    ! to Bohr, as done in the prep module.
    ! 
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  27 Oct 2011    A. Reid        Development of subroutine
    !  24 May 2012    A. Reid        Parallelization with MPI
    !
    ! =======================================================================

    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================

    ! Calling variables:
    double precision,DIMENSION(d,p,n),INTENT(IN) :: xin	! Ring polymer geometry
    double precision,DIMENSION(d*p*n) :: x			! Ring polymer geometry

    ! Local variables:
    INTEGER :: xdim					! Dimension of xbead
    double precision,DIMENSION(:),ALLOCATABLE :: xbead	! Bead coordinates, to use in calculating potential
    double precision :: beadpot				! To receive calculated potential
    INTEGER :: i					! Index of time slices
    INTEGER :: ind					! Index of DO loops

    ! MPI variables:
    double precision,DIMENSION(d*p,n) :: beadsendbuf	! Linear polymer geometry split into n beads, to be scattered to processes
    double precision,DIMENSION(:),ALLOCATABLE :: beadrecvbuf	! Bead geometries, to be received by one of nprocs processes
    double precision,DIMENSION(d*p) :: beadgeom		! Bead geometry of a single bead for calculation of potential
    double precision,INTENT(OUT) :: potsum			! Sum of single bead potentials, in the root process
    double precision :: potsendbuf                        	! Single bead potential, calculated in one of nprocs processes
    double precision,DIMENSION(n) :: potrecvbuf		! Array of pots for all beads, to be gathered from nprocs processes
    double precision :: potential				! Sum of single bead potentials, non-zero except in the root process
    INTEGER :: error, nprocs, namelen, myrank

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
        CALL threetoone ( d, p, n, xin, x )
        DO i = 1, n
            beadsendbuf(1:(d*p),i) = x( ((i-1)*(d*p))+1 : (i*(d*p)) )
        END DO
    END IF initialsetupmulti

    ALLOCATE(beadrecvbuf(d*p*n/nprocs))

    ! Initial setup of the x array is complete, so we now scatter the beads to the different processes
    CALL MPI_Scatter(beadsendbuf, d*p*n/nprocs, MPI_DOUBLE_PRECISION, &
         beadrecvbuf, d*p*n/nprocs, MPI_DOUBLE_PRECISION, root, &
         MPI_COMM_WORLD, error)

    ! The code here will be executed in each process
    ! Setting terms to zero
    potsendbuf = 0.0d0
    beadpot = 0.0d0

    DO i = 1, n/nprocs
        beadgeom = beadrecvbuf( ((i-1)*(d*p))+1 : (i*(d*p)) )
        call pes_obj%pot(beadgeom,beadpot)
        beadpot=beadpot-(cdata%potential)                       ! Rezeroing of potential (relative to minimum potential)
        potsendbuf = potsendbuf + beadpot			! Running total of potential for all beads in the ring polymer
    END DO

    ! Calculation of bead potentials is complete, so we now sum/broadcast them to all processes
    CALL MPI_Allreduce(potsendbuf, potsum, 1, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_COMM_WORLD, error)

    END SUBROUTINE extpot

! ===========================================================================================================

    SUBROUTINE totpotx ( xin, totalpot )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To calculate the total ring-polymer potential from using subroutines
    ! for the internal and external parts.  Takes input in format compatible
    ! with Marko's L-BFGS code.
    !
    ! 
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  21 Nov 2011    A. Reid        Development of subroutine
    !
    ! =======================================================================

    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================

    ! Calling variables:
    double precision,DIMENSION(d*p*n),INTENT(IN) :: xin               ! Linear polymer geometry in 1D array
    double precision,INTENT(OUT) :: totalpot                          ! Total ring polymer potential

    ! Local variables:
    double precision,DIMENSION(d,p,n) :: x                  ! Linear polymer geometry
    double precision,DIMENSION(d,p,n) :: xscaled             ! Scaled linear polymer geometry
    double precision :: externalpot                          ! External ring polymer potential
    double precision :: internalpot                          ! Internal ring polymer potential

    ! =======================================================================
    ! CONVERTING FROM 1D ARRAY TO 3D ARRAY
    ! =======================================================================

    CALL onetothree ( d, p, n, xin, x )

    ! =======================================================================
    ! RP POTENTIAL ('INTERNAL' PART)
    ! =======================================================================

    CALL intpotscaled ( x, internalpot )

    ! =======================================================================
    ! RP POTENTIAL ('EXTERNAL' PART)
    ! =======================================================================

    CALL extpot ( x, externalpot )

    ! =======================================================================
    ! RP POTENTIAL (TOTAL)
    ! =======================================================================

    totalpot = internalpot + externalpot

    END SUBROUTINE totpotx

! ===========================================================================================================

!      SUBROUTINE onebeadpotx ( xonebead, onebeadpotential )
!      !
!      ! =======================================================================
!      ! PURPOSE OF THIS SUBROUTINE
!      ! =======================================================================
!      ! To calculate the 'external' potential for a single bead.  Note that the
!      ! input geometry must have already been converted from Angstroms to Bohr,
!      ! as done in the prep module.
!      !
!      ! Record of revisions:
!      !
!      !  Date           Person         Change implemented
!      !  ----           ------         ------------------
!      !  5 Dec 2011     A. Reid        Development of subroutine
!      !
!      ! =======================================================================
!  
!      ! =======================================================================
!      ! DECLARATION OF VARIABLES
!      ! =======================================================================
!  
!      ! Calling variables:
!      REAL(kind=dp),DIMENSION(d*p),INTENT(IN) :: xonebead                 ! Ring polymer geometry
!      REAL(kind=dp),INTENT(OUT) :: onebeadpotential                          ! External ring polymer potential
!  
!      ! Local variables:
!      REAL(kind=dp),DIMENSION(d*p) :: xcalc                  
!  
!      xcalc=xonebead
!      onebeadpotential=f(xcalc)                        ! See pes_shell for function
!  
!      END SUBROUTINE onebeadpotx

! ===========================================================================================================

    SUBROUTINE pathway(x,pathwayfile)
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To produce a file called "potential" which can make a plot of V(x) vs
    ! distance along the instanton pathway.
    ! 
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  9 Jul 2012     J. Richardson  Development of subroutine
    !  9 Jul 2012     A. Reid        Tweaks to filename and layout
    !
    ! =======================================================================

    double precision,DIMENSION(d,p,n),INTENT(IN) :: x	! Linear polymer geometry in 3D array
    double precision,dimension(d,p,n) :: xpes
    double precision,dimension(p) :: mpes
    CHARACTER(len=64),INTENT(IN) :: pathwayfile		! File to save instanton pathway
    double precision :: V, r ! potential and mass-weighted distance along pathway
    integer :: i

    xpes = x
    mpes = masses
    if (.not. transformed_to_pot) then
            !transform first bead and masses to pes coordinates:
            call pes_obj%vec_to_pes(xpes(:,:,1))
            call pes_obj%atomvec_to_pes(mpes)
    end if


    r = 0.0d0						! Set mass-weighted distance to zero at start
    open(21,file=pathwayfile)
    call pes_obj%pot(xpes(:,:,1), V)                       ! Calculate potential
    write(21,*) r, V					! Write pair of data points to file
    do i = 2, n
        !transform bead to pes coordinates:
        if (.not. transformed_to_pot) call pes_obj%vec_to_pes(xpes(:,:,i))
        call pes_obj%pot(xpes(:,:,i), V)                   ! Calculate potential
        r = r + sqrt(sum(mpes(:)*sum((xpes(:,:,i)-xpes(:,:,i-1))**2,1)))	! Increment mass-weighted distance
        write(21,*) r, V				! Write pair of data points to file
    end do
    close(21)

    END SUBROUTINE pathway

! ===========================================================================================================

END MODULE pots
