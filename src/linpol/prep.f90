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
MODULE prep
USE parameters
USE phys_const
!
! =======================================================================
! PURPOSE OF THIS MODULE
! =======================================================================
! To contain the subroutines to prepare for a linear-polymer minimization
! and/or calculation of a hessian matrix and accompanying eigenvalues.
! 
! Contents:
!     xyzread
!     initialends
!     singlebead
!     array
!     data
!     externalevals
!
! Record of revisions:
!
!  Date           Person         Change implemented
!  ----           ------         ------------------
!  22 Nov 2011    A. Reid        Development of module
!  29 Nov 2011    A. Reid        Conversion to generate variables to share
!  9 May 2012     A. Reid        Redundant subroutines deleted
!  4 Jul 2013     A. Reid        Added externalevals subroutine
!
! =======================================================================

IMPLICIT NONE
SAVE                                       ! Ensures all declared data values will be preserved

! =======================================================================
! DECLARATION OF GLOBAL VARIABLES AVAILABLE TO OTHER PROGRAM UNITS
! =======================================================================

! Values fixed for all uses of the code:
INTEGER :: d = 3                           ! No. of dimensions for coding purposes

! Values to be generated from the code in this module:
INTEGER :: actualdim                       ! No. of dimensions in actual system
INTEGER :: ninit                           ! No. of time slices / f-dimensional beads as calculated from initialisation subroutine
INTEGER :: n                               ! No. of time slices / f-dimensional beads
INTEGER :: p                               ! No. of atoms
double precision :: beta                      ! Beta = 1/kT
double precision :: betan                     ! Betan = beta/n
INTEGER :: dof                             ! No. of degrees of freedom
CHARACTER(LEN=256) :: clustersfile			! File containing the minimum cluster potentials

! Arrays to be generated from the code in this module:

double precision,DIMENSION(:),ALLOCATABLE :: masses            ! Masses
CHARACTER(len=2),DIMENSION(:),ALLOCATABLE :: atomtags       ! Atom tags
double precision,DIMENSION(:,:),ALLOCATABLE :: xsingle         ! Starting geometry
logical :: transformed_to_pot = .false.                     ! have the coordinates been transformed to convention used by the potential?
double precision,DIMENSION(:,:,:),ALLOCATABLE :: x             ! Linear-polymer
double precision,DIMENSION(:,:),ALLOCATABLE :: xA, xB          ! fixed ends
!REAL(kind=dp) :: clusterminpots(100)		! Array to contain the minimum potentials, with index = clustersize
!REAL(kind=dp) :: zeropot
double precision,DIMENSION(:),ALLOCATABLE :: extevals             ! For debugging/testing purposes

! =======================================================================
! DECLARATION OF LOCAL VARIABLES
! =======================================================================

INTEGER :: ind                                        ! Index of DO loops
INTEGER :: i                                          ! Index of time slices


CONTAINS

! ===========================================================================================================

    SUBROUTINE xyzread ( filename )
    USE MPI
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To read an xyz file and hence determine the number of atoms, plus beta.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  22 Dec 2011    A. Reid        Development of subroutine
    !  20 Nov 2012    M. Herbst      Minor changes to incorporate commandline option for beta
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================

    ! Calling variables:
    CHARACTER(len=64),INTENT(IN) :: filename		! File to open
    
    ! Local variables:
    CHARACTER(80) :: msg				! Status message
    INTEGER :: status					! IO status message
    INTEGER :: myrank, error				! For MPI

    ! =======================================================================
    ! GENERATE MPI DATA
    ! =======================================================================
    CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)

    ! =======================================================================
    ! PROCEED WITH SUBROUTINE
    ! =======================================================================

    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    WRITE(99,'(A)') 'Reading atomic coordinates from xyz file'
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    
    
    OPEN (UNIT=1, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=status, IOMSG=msg )
    ! Note that writing UNIT='1' (with inverted commas) would cause the program to fail!
    openif: IF ( status == 0 ) THEN
        !OPEN succeeded - proceed with program
    
        READ ( 1, *, IOSTAT=status ) p			! Read number of atoms from first line in file
        pread: IF ( status /= 0 ) THEN
            IF (myrank == root) THEN			! To avoid duplication in non-root processes
                WRITE(99,*) '********************'				! Error check
                WRITE(*,*) '********************'				! Error check
                WRITE(99,*) 'Error reading p!'					! Error check
                WRITE(*,*) 'Error reading p!'					! Error check
                WRITE(99,*) 'Check first line of xyz file used as input'	! Error check
                WRITE(*,*) 'Check first line of xyz file used as input'		! Error check
                WRITE(99,*) '********************'				! Error check
                WRITE(*,*) '********************'				! Error check
            END IF					! To avoid duplication in non-root processes
            STOP
        ELSE
            IF (myrank == root) THEN			! To avoid duplication in non-root processes
                WRITE(99,1040) p						! Echo p to user
                WRITE(*,1040) p						! Echo p to user
                1040 FORMAT ('No. of atoms (p) is ', I3)			! Echo p to user
            END IF
        ENDIF pread

        IF (myrank == root) THEN			! To avoid duplication in non-root processes
            cluster: IF ( p == 6 ) THEN
                WRITE(99, '(A)' ) 'This is a water dimer'
                WRITE(*, '(A)' ) 'This is a water dimer'                                                                            
            ELSE IF ( p == 9 ) THEN
                WRITE(99, '(A)' ) 'This is a water trimer'
                WRITE(*, '(A)' ) 'This is a water trimer'    
            ELSE IF ( p == 12 ) THEN
                WRITE(99, '(A)' ) 'This is a water tetramer'
                WRITE(*, '(A)' ) 'This is a water tetramer'
            ELSE IF ( p == 15 ) THEN
                WRITE(99, '(A)' ) 'This is a water pentamer'
                WRITE(*, '(A)' ) 'This is a water pentamer'
            ELSE IF ( p == 18 ) THEN
                WRITE(99,'(A)') 'This is a water hexamer'
                WRITE(*,'(A)') 'This is a water hexamer'
            ELSE IF ( p == 21 ) THEN
                WRITE(99,'(A)') 'This is a water heptamer'
                WRITE(*,'(A)') 'This is a water heptamer'
            ELSE IF ( p == 24 ) THEN
                WRITE(99,'(A)') 'This is a water octamer'
                WRITE(*,'(A)') 'This is a water octamer'
            ELSE IF ( p == 27 ) THEN
                WRITE(99,'(A)') 'This is a water nonamer'
                WRITE(*,'(A)') 'This is a water nonamer'
            ELSE IF ( p == 30 ) THEN
                WRITE(99,'(A)') 'This is a water decamer'
                WRITE(*,'(A)') 'This is a water decamer'
            ELSE
                WRITE(99,'(A,I2,A)') 'This is a water ', p/3, '-mer'
                WRITE(*,'(A,I2,A)') 'This is a water ', p/3, '-mer'
            ENDIF cluster
        ENDIF						! To avoid duplication in non-root processes
    

        !edit by mfh for commandline change of beta:
        READ ( 1, *, IOSTAT=status ) beta					! Read beta from second line in file
        if ( masterbeta .gt. 0.0 ) then 
                !note previous line is done in either case to make sure program
                !is in consistent (can't be bothered to check if the current
                !line actually needs to be pushed forward for the rest of the
                !program.
                beta = masterbeta
                IF (myrank == root) THEN			! To avoid duplication in non-root processes
                     WRITE(99,1050) beta						! Echo beta to user
                     WRITE(*,1050) beta						! Echo beta to user
                ENDIF					! To avoid duplication in non-root processes
        else
                betaread: IF ( status /= 0 ) THEN
                    IF (myrank == root) THEN			! To avoid duplication in non-root processes
                        WRITE(99,*) '********************'				! Error check
                        WRITE(*,*) '********************'				! Error check
                        WRITE(99,*) 'Error reading beta!'				! Error check
                        WRITE(*,*) 'Error reading beta!'				! Error check
                        WRITE(99,*) 'Check second line of xyz file used as input'	! Error check
                        WRITE(*,*) 'Check second line of xyz file used as input'	! Error check
                        WRITE(99,*) '********************'				! Error check
                        WRITE(*,*) '********************'				! Error check
                    ENDIF					! To avoid duplication in non-root processes
                    STOP
                ELSE
                    IF (myrank == root) THEN			! To avoid duplication in non-root processes
                        WRITE(99,1050) beta						! Echo beta to user
                        WRITE(*,1050) beta						! Echo beta to user
                        1050 FORMAT ('beta is ', F12.3)					! Echo beta to user
                    ENDIF					! To avoid duplication in non-root processes
                ENDIF betaread
        end if 


    ELSE openif
        !OPEN failed
        WRITE(99,1000) status
        1000 FORMAT ('Error opening file: IOSTAT = ', I6)
        WRITE(99,*) TRIM(msg)
        WRITE(99,*) 'FATAL ERROR!'
        ! Repeated for printing to screen:
        WRITE(*,1000) status
        WRITE(*,*) TRIM(msg)
        WRITE(*,*) 'FATAL ERROR!'
        STOP
    ENDIF openif


    END SUBROUTINE xyzread

! ===========================================================================================================

    SUBROUTINE initialends ( filename )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To read an xyz file and hence determine the orientation of the system
    ! being studied.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  27 Oct 2011    A. Reid        Development of 'initial' subroutine
    !  2 Nov 2011     A. Reid        Adaptation into 'initialends' subroutine
    !  22 Dec 2011    A. Reid        File read functionality moved to 'xyzread'
    !  22 Dec 2011    A. Reid        Merged with old 'initial' subroutine
    !  1 June 2012    A. Reid        Turned off printing of beads to log file
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Calling variables:
    CHARACTER(len=64),INTENT(IN) :: filename              ! File to open
    
    ! Local variables:
    CHARACTER(80) :: msg                                  ! Status message
    INTEGER :: status                                     ! IO status message
    INTEGER :: nvals = 0                                  ! No. of values read in from file
    CHARACTER(2) :: atom                                  ! Value read in from file
    double precision :: xcoord                               ! Value read in from file
    double precision :: ycoord                               ! Value read in from file
    double precision :: zcoord                               ! Value read in from file
    
    
    CALL xyzread ( filename )

    ! =======================================================================
    ! TESTING NUMBER OF VARIABLES SUPPLIED BY USER
    ! =======================================================================
    
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    WRITE(99,'(A)') 'Testing number of dimensions'
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    
                READ ( 1, *, IOSTAT=status ) atom, xcoord, ycoord, zcoord                  ! Read first line of values
    
                IF ( ( zcoord == 0 ).AND.( ycoord == 0 ).AND.( xcoord == 0 ).AND.( status /= 0 ) )  THEN
                    WRITE(99,*) 'x, y and z co-ordinates not found'
                    WRITE(99,*) 'FATAL ERROR!'
                    STOP
                    actualdim = 0
                ELSEIF ( ( zcoord == 0 ).AND.( ycoord == 0 ).AND.( status /= 0 ) ) THEN
                    WRITE(99,*) 'y and z co-ordinates not found'
                    WRITE(99,*) 'One-dimensional system'
                    actualdim = 1
                ELSEIF ( ( zcoord == 0 ).AND.( status /= 0 ) ) THEN
                    WRITE(99,*) 'z co-ordinate not found'
                    WRITE(99,*) 'Two-dimensional system'
                    actualdim = 2
                ELSE
                    WRITE(99,*) 'x, y and z co-ordinates found'
                    WRITE(99,*) 'Three-dimensional system'
                    actualdim = 3
                ENDIF
    
    REWIND (1, IOSTAT=status)                                        ! Rewind the xyz file to start from scratch
    READ ( 1, *, IOSTAT=status )                                                       ! Read line but do not record
    READ ( 1, *, IOSTAT=status )                                                       ! Read line but do not record

    ! =======================================================================
    ! READING COORDINATE VALUES FROM XYZ FILE
    ! =======================================================================
    
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    WRITE(99,'(A)') 'Reading values from file'
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    
    dimsread: IF (actualdim == 3) THEN
    ! ----------------------
    ! THREE-DIMENSIONAL CASE
    ! ----------------------

        i = 1
    
        threedreadloop: DO
    
            !WRITE(99,1060) i
            1060 FORMAT ('f-dimensional bead number: ' I6)
    
            threedreadsubloop: DO ind = 1, p, 1                                              ! DO loop from 1 to p (no. of atoms) in increments of 1
                READ ( 1, *, IOSTAT=status ) atom, xcoord, ycoord, zcoord                  ! Read next line of values
                IF ( status /= 0 ) EXIT                                                    ! Negative means end reached, positive (non-zero) means error
                nvals = nvals + 1                                                          ! Valid, hence increase counter of no. of values
                !WRITE(99,1010) atom, xcoord, ycoord, zcoord            ! Echo line number and values to screen
                1010 FORMAT (A2, ' ', F7.3, ' ', F7.3, ' ', F7.3)                             ! Echo line number and values to screen
            END DO threedreadsubloop
    
            READ ( 1, *, IOSTAT=status )                                                       ! Read line but do not record
            IF ( status /= 0 ) EXIT                                                            ! Negative means end reached, positive (non-zero) means error
            READ ( 1, *, IOSTAT=status )                                                       ! Read line but do not record
            IF ( status /= 0 ) EXIT                                                            ! Negative means end reached, positive (non-zero) means error
    
            i = i + 1                                                               ! There are more atoms, hence increase counter of no. of beads
    
        END DO threedreadloop
    
        ninit = i                              ! Now loop has completed, assign counted value of i to variable ninit
    
    
    ELSEIF (actualdim == 2) THEN
    ! ----------------------
    ! TWO-DIMENSIONAL CASE
    ! ----------------------    

        i = 1
    
        twodreadloop: DO
    
            WRITE(99,1160) i
            1160 FORMAT ('f-dimensional bead number: ' I6)
    
            twodreadsubloop: DO ind = 1, p, 1                                              ! DO loop from 1 to p (no. of atoms) in increments of 1
                READ ( 1, *, IOSTAT=status ) atom, xcoord, ycoord                  ! Read next line of values
                zcoord = 0.0
                IF ( status /= 0 ) EXIT                                                    ! Negative means end reached, positive (non-zero) means error
                nvals = nvals + 1                                                          ! Valid, hence increase counter of no. of values
                !WRITE(99,1110) atom, xcoord, ycoord, zcoord            ! Echo line number and values to screen
                1110 FORMAT (A2, ' ', F7.3, ' ', F7.3, ' ', F7.3)                             ! Echo line number and values to screen
            END DO twodreadsubloop
    
            READ ( 1, *, IOSTAT=status )                                                       ! Read line but do not record
            IF ( status /= 0 ) EXIT                                                            ! Negative means end reached, positive (non-zero) means error
            READ ( 1, *, IOSTAT=status )                                                       ! Read line but do not record
            IF ( status /= 0 ) EXIT                                                            ! Negative means end reached, positive (non-zero) means error
    
            i = i + 1                                                               ! There are more atoms, hence increase counter of no. of beads
    
        END DO twodreadloop
    
        ninit = i                              ! Now loop has completed, assign counted value of i to variable ninit
    
    
    ELSEIF (actualdim == 1) THEN
    ! ----------------------
    ! ONE-DIMENSIONAL CASE
    ! ----------------------
    
        i = 1
    
        onedreadloop: DO
    
            WRITE(99,1260) i
            1260 FORMAT ('f-dimensional bead number: ' I6)
    
            onedreadsubloop: DO ind = 1, p, 1                                              ! DO loop from 1 to p (no. of atoms) in increments of 1
                READ ( 1, *, IOSTAT=status ) atom, xcoord                  ! Read next line of values
                zcoord = 0.0
                ycoord = 0.0
                IF ( status /= 0 ) EXIT                                                    ! Negative means end reached, positive (non-zero) means error
                nvals = nvals + 1                                                          ! Valid, hence increase counter of no. of values
                !WRITE(99,1210) atom, xcoord, ycoord, zcoord            ! Echo line number and values to screen
                1210 FORMAT (A2, ' ', F7.3, ' ', F7.3, ' ', F7.3)                             ! Echo line number and values to screen
            END DO onedreadsubloop
    
            READ ( 1, *, IOSTAT=status )                                                       ! Read line but do not record
            IF ( status /= 0 ) EXIT                                                            ! Negative means end reached, positive (non-zero) means error
            READ ( 1, *, IOSTAT=status )                                                       ! Read line but do not record
            IF ( status /= 0 ) EXIT                                                            ! Negative means end reached, positive (non-zero) means error
    
            i = i + 1                                                               ! There are more atoms, hence increase counter of no. of beads
    
        END DO onedreadloop
    
        ninit = i                              ! Now loop has completed, assign counted value of i to variable ninit
    
    ENDIF dimsread
    
    ! ----------------------
    ! ERROR CHECKING
    ! ----------------------
    ! The DO loop has completed, which could be because status did not equal zero at some point,
    ! i.e. the EXIT command was triggered by a READ error, or it could be because the end of the
    ! file was reached.  We now investigate:
    
        readif: IF ( status > 0 ) THEN                ! Read error occurred
            WRITE(99,1020) nvals + 1
            1020 FORMAT ('An error occurred reading line ', I6)
        ELSE readif                                   ! End of data reached
            WRITE(99,*)                                    ! Create space in output
            WRITE(99,1030) nvals
            1030 FORMAT ('End of file reached. There were ', I6, ' atom-containing lines in the file.')
            WRITE(99,1070) ninit
            1070 FORMAT ('These lines were split into ', I6, ' f-dimensional beads.')
    
            check: IF ( nvals / ninit /= p ) THEN                    ! Error check
                    WRITE(99,'(A)')
                    WRITE(99,'(A)') '*******************************************************************'
                    WRITE(99,'(A)') 'ERROR: (no. of lines) / (no. of beads) does NOT equal no. of atoms!'
                    WRITE(99,'(A)') '*******************************************************************'
                    WRITE(99,'(A)')
                ELSE
                    WRITE(99,'(A)') 'CHECK: (No. of lines) / (no. of beads) equals no. of atoms, as expected.'                                                    
            ENDIF check
    
    
        ENDIF readif
    
    
    CLOSE (UNIT=1)

    END SUBROUTINE initialends

! ===========================================================================================================

    SUBROUTINE singlebead ( filename )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To read an xyz file and hence determine the parameters of the system
    ! being studied.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  4 Nov 2011     A. Reid        Adaptation from 'initialends' subroutine
    !  18 Mar 2012    A. Reid        Added capability to deal with a single
    !                                bead, and edited out old comments.
    !                                Original subroutine saved as startendold
    !  9 May 2012     A. Reid        Edited out old comments
    !  9 May 2012     A. Reid        Got rid of scaling procedures and 1D/2D
    !  9 May 2012     A. Reid        Name changed from 'startend' to 'singlebead'
    !
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Calling variables:
    CHARACTER(len=64),INTENT(IN) :: filename		! File to open
    
    ! Local variables:
    CHARACTER(80) :: msg				! Status message
    INTEGER :: status					! IO status message
    CHARACTER(4) :: atomtag				! Value read in from file
    double precision :: xcoord				! Value read in from file
    double precision :: ycoord				! Value read in from file
    double precision :: zcoord				! Value read in from file
    INTEGER :: atom					! Index of time slices
    double precision, DIMENSION(d,p) :: xsinglegen		! Starting geometry
    double precision, DIMENSION(p) :: massesgen		! Masses
    CHARACTER(len=4),DIMENSION(p) :: atomtagsgen	! Atom tags

    
    
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    WRITE(99,'(A)') 'Reading in single-bead data from xyz file'
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    
    OPEN (UNIT=1, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=status, IOMSG=msg )
    
    READ ( 1, *, IOSTAT=status )			! Read line but do not record
    READ ( 1, *, IOSTAT=status )			! Read line but do not record
    
    i = 1
    
    readloop: DO
    
        readsubloop: DO atom = 1, p, 1				! DO loop from 1 to p (no. of atoms) in increments of 1
            READ ( 1, *, IOSTAT=status ) atomtag, xcoord, ycoord, zcoord	! Read next line of values
            IF ( status /= 0 ) EXIT						! Negative means end reached, positive (non-zero) means error
    
            IF ( atomtag(1:1) == "H" ) massesgen( atom ) = auamu*Hamu
            IF ( atomtag(1:1) == "D" ) massesgen( atom ) = auamu*Damu
            IF ( atomtag(1:1) == "O" ) massesgen( atom ) = auamu*Oamu

            atomtagsgen(atom) = atomtag(1:4)
    
            xsinglegen(1,atom) = xcoord
            xsinglegen(2,atom) = ycoord
            xsinglegen(3,atom) = zcoord
    
        END DO readsubloop
    
        READ ( 1, *, IOSTAT=status )					! Read line but do not record
        IF ( status /= 0 ) EXIT						! Negative means end reached, positive (non-zero) means error
        READ ( 1, *, IOSTAT=status )					! Read line but do not record
        IF ( status /= 0 ) EXIT						! Negative means end reached, positive (non-zero) means error
    
    END DO readloop
    

    CLOSE (UNIT=1)


    ALLOCATE(xsingle(d,p))				! Allocate shared array to hold starting geometry
    ALLOCATE(masses(p))					! Allocate shared array to hold masses
    ALLOCATE(atomtags(p))				! Allocate shared array to hold atom tags

    xsingle = xsinglegen / auang				! Copy calculated starting geometry to shared array and convert from Angstroms to bohr
    masses = massesgen					! Copy generated masses to shared array
    atomtags = atomtagsgen				! Copy generated atom tags to shared array

    END SUBROUTINE singlebead

! ===========================================================================================================

    SUBROUTINE array ( filename )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! Populates the 3-D array and masses array based on the parameters
    ! provided.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  27 Oct 2011    A. Reid        Development of subroutine
    !
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Calling variables:
    CHARACTER(len=64),INTENT(IN) :: filename              ! File to open   

    ! Local variables:
!    CHARACTER(2) :: atom                                   ! Value read in from file
    double precision :: xcoord                                         ! Value read in from file
    double precision :: ycoord                                         ! Value read in from file
    double precision :: zcoord                                         ! Value read in from file
!    INTEGER :: ind                                         ! Index of DO loops
    INTEGER :: i                                           ! Index of time slices
    INTEGER :: status                                     ! IO status message
!added:
    double precision, DIMENSION(p) :: massesgen             ! Masses
    CHARACTER(len=4),DIMENSION(p) :: atomtagsgen      ! Atom tags
    CHARACTER(4) :: atomtag                                  ! Value read in from file
    INTEGER :: atom                                   ! Index of atoms
    CHARACTER(80) :: msg                                  ! Status message
    double precision, DIMENSION(:,:,:),ALLOCATABLE :: xgen             ! Linear-polymer

    ALLOCATE(xgen(d,p,n))


    OPEN (UNIT=1, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=status, IOMSG=msg )
    ! Note that writing UNIT='1' (with inverted commas) would cause the program to fail!
    openif: IF ( status == 0 ) THEN
        !OPEN succeeded - proceed with program
    !WRITE(99,*) 'DEBUG: OPEN SUCCESSFUL'
        
    ELSE openif
        !OPEN failed
        WRITE(99,1000) status
        1000 FORMAT ('Error opening file: IOSTAT = ', I6)
        WRITE(99,*) TRIM(msg)
    ENDIF openif


    
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    WRITE(99,'(A)') "Populating the 3-D array"
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    
    REWIND (1, IOSTAT=status)                                        ! Rewind the xyz file to start from scratch
    i = 1                                                            ! Reset counter to 1
    
    readloop2: DO
    
        READ ( 1, *, IOSTAT=status )                                                       ! Read line but do not record
        IF ( status /= 0 ) EXIT                                                            ! Negative means end reached, positive (non-zero) means error
        READ ( 1, *, IOSTAT=status )                                                       ! Read line but do not record
        IF ( status /= 0 ) EXIT                                                            ! Negative means end reached, positive (non-zero) means error
    
        readsubloop2: DO atom = 1, p, 1                                              ! DO loop from 1 to p (no. of atoms) in increments of 1
            READ ( 1, *, IOSTAT=status ) atomtag, xcoord, ycoord, zcoord                  ! Read next line of values
            IF ( status /= 0 ) EXIT                                                    ! Negative means end reached, positive (non-zero) means error
            xgen( 1, atom, i ) = xcoord
            xgen( 2, atom, i ) = ycoord
            xgen( 3, atom, i ) = zcoord
    
            IF ( atomtag(1:1) == "H" ) massesgen( atom ) = auamu*Hamu
            IF ( atomtag(1:1) == "D" ) massesgen( atom ) = auamu*Damu
            IF ( atomtag(1:1) == "O" ) massesgen( atom ) = auamu*Oamu

            atomtagsgen(atom) = atomtag(1:4)
    
        END DO readsubloop2
    
            i = i + 1                                                               ! There are more atoms, hence increase counter of no. of beads
    
    END DO readloop2
    
    ALLOCATE( x( d, p, n ) )	! Create the blank 3-D array
    x = xgen / auang		! Copy calculated starting geometry to shared array and convert from Angstroms to bohr.
				! Lengths now in bohr for rest of the program, except when they need to be converted into
				! the native units of the potential being utilised in order to make a potential call; this
				! is dealt with by pes_interface.f90
   
    ALLOCATE(masses(p))        ! Allocate shared array to hold masses
    masses = massesgen         ! Copy generated masses to shared array

    ALLOCATE(atomtags(p))      ! Allocate shared array to hold atom tags
    atomtags = atomtagsgen     ! Copy generated atom tags to shared array


    END SUBROUTINE array

! ===========================================================================================================

    SUBROUTINE data ( )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! Reports data on the 3-D array and calculates betan.
    !
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  27 Oct 2011    A. Reid        Development of subroutine
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Local variables:
    INTEGER :: entries                         ! No. of entries in 3-D array, x
   
    
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    WRITE(99,'(A)') "Data for the 3-D array"
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    
    dof = d * p                                  ! Number of degrees of freedom

    SELECT CASE(endchoice)
        CASE(1)  ! Fixed ends
        !    betan = beta / (n+1)
            betan = beta / n
        CASE(2)  ! Unfixed ends
        !    betan = beta / (n-1)
            betan = beta / n
        CASE DEFAULT  ! endchoice not specified - defaults to unfixed ends
        !    betan = beta / (n-1)
            betan = beta / n
    END SELECT
    
    entries=d*p*n                                ! Calculate number of entries in 3-D array
    
    WRITE(99, 1080 ) dof                        ! Write number of degrees of freedom
    1080 FORMAT ('f = ', I6)
    WRITE(99, 1090 ) betan                      ! Write betaN
    1090 FORMAT ('betan = ', F12.3)
    
    
    WRITE(99,1100) entries                        ! Write number of entries in 3-D array (x)
    1100 FORMAT ('d*p*n = no. of entries in 3-D array (called x) =', I6)
    
    
    END SUBROUTINE data

! ===========================================================================================================

    SUBROUTINE externalevals ( )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To read eigenvalues from a file for debugging/testing purposes.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  4 Jul 2013     A. Reid        Development of subroutine
    !
    ! =======================================================================

    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Local variables:
    CHARACTER(len=512) :: filename              ! File to open
    double precision :: firsteval                                         ! Value read in from file
    double precision :: secondeval                                         ! Value read in from file
    double precision :: thirdeval                                         ! Value read in from file
    INTEGER :: i                                           ! Index of lines in file
    INTEGER :: status                                     ! IO status message
    CHARACTER(80) :: msg                                  ! Status message
    double precision, DIMENSION(:),ALLOCATABLE :: evalsgen             ! eigenvalues

    filename = 'eigvals'
    ALLOCATE(evalsgen(d*p*n))

    OPEN (UNIT=88, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=status, IOMSG=msg )
    openif: IF ( status == 0 ) THEN
        !OPEN succeeded - proceed with program
    !WRITE(99,*) 'DEBUG: OPEN SUCCESSFUL'
        
    ELSE openif
        !OPEN failed
        WRITE(99,1000) status
        1000 FORMAT ('Error opening file: IOSTAT = ', I6)
        WRITE(99,*) TRIM(msg)
    ENDIF openif

   
    REWIND (88, IOSTAT=status)                                        ! Rewind the xyz file to start from scratch
    
        readsubloop2: DO i = 1, d*p*n/3                                              
            READ ( 88, *, IOSTAT=status ) firsteval, secondeval, thirdeval                  ! Read next line of values
            IF ( status /= 0 ) EXIT                                                    ! Negative means end reached, positive (non-zero) means error
            evalsgen( ((i-1)*3) + 1 ) = firsteval
            evalsgen( ((i-1)*3) + 2 ) = secondeval
            evalsgen( ((i-1)*3) + 3 ) = thirdeval
   
        END DO readsubloop2
        
    ALLOCATE( extevals( d*p*n ) )					! Create the blank extevals array
    extevals = evalsgen			! Copy recorded evals to shared array
   
    END SUBROUTINE externalevals

! ===========================================================================================================

END MODULE prep
