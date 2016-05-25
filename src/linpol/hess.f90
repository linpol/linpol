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
MODULE hess
USE parameters						! Module containing parameters to be pulled into all other modules, e.g. commenting information
USE prep						! Linear polymer preparation module
IMPLICIT NONE
!
! =======================================================================
! PURPOSE OF THIS MODULE
! =======================================================================
! To contain the subroutines to calculate the hessian matrix of second
! derivatives, to diagonalize it and to output it to a file.
! 
! Contents:
!     hesscalc
!     hesscalcsingle
!     hessext
!     hessout
!     diag
!     diagband
!     diagsymm
!
! Record of revisions:
!
!  Date           Person         Change implemented
!  ----           ------         ------------------
!  22 Nov 2011    A. Reid        Development of module
!  14 Feb 2012    A. Reid        Changed finite diffs to asymmetrical
!  14 Feb 2012    A. Reid        Changed Hessian to banded storage
!  9 May 2012     A. Reid        Removed old commented-out code
!  9 May 2012     A. Reid        Format descriptors tidied up and IF comments removed
!  24 Nov 2012    M. Herbst      Parellelised Hessian calculation
!
! =======================================================================

CONTAINS

! ===========================================================================================================

    SUBROUTINE hesscalc ( xmin, hessian, useMPI)
    USE funcs						! Miscellaneous functions module
    USE grads						! Gradient calculation module
    USE mpi
    use parameters, only: root
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To take an an xyz file representing a linear polymer and calculate the
    ! Hessian matrix of second derivatives.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  7 Nov 2011     A. Reid        Development of subroutine
    !  15 Nov 2011    A. Reid        Completion of initial version
    !  20 Jan 2012    A. Reid        Updated: shared variables / non-scaled x
    !  31 Jan 2012    A. Reid        Added debugging code
    !  27 Feb 2012    A. Reid        Grad calculations moved to subroutines
    !  31 Mar 2012    J. Richardson  Simplified
    !  24 Nov 2012    M. Herbst      Parallelised hessian calculation
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Calling variables:
    double precision,DIMENSION(d,p,n),INTENT(IN) :: xmin			! Minimum geometry
    logical,intent(in) :: useMPI
    double precision,DIMENSION(((d*p)+1),d*p*n),INTENT(INOUT) :: hessian	! fixed-ended mass-weighted linear polymer hessian (upper) (band storage)
    
    ! Local variables (user may wish to alter the assigned values):
    INTEGER :: terminalint = 10000					! Interval between writes to terminal
    double precision :: delta = 0.0001					! Finite difference in co-ordinate

    ! Local variables (no assigned values):
    INTEGER :: beadi							! Index of time slices
    INTEGER :: ib, j, k							! Indices to loop over to create banded Hessian
    INTEGER :: totalentries						! Total no. of entries in Hessian
    INTEGER :: entriestocalc						! Total no. of entries to be calculated
    INTEGER :: entrycount						! Hessian entry no.
    double precision :: percentcomplete					! Percentage completion
    double precision :: totalentriesreal				! Total no. of entries in Hessian (as a real number)
    double precision :: entrycountreal					! Hessian entry no. (as a real number)
    REAL :: terminalreal						! For entrycount test for writing to terminal
    INTEGER :: fileint							! For entrycount test for writing to file
    REAL :: filereal							! For entrycount test for writing to file
    double precision,DIMENSION(d*p,d*p) :: H				! Hessian of one bead

    character*64 :: hessfile  !-------------------------------------------------------------------------------------------------------------------

    !MPI setup
    INTEGER :: error, nprocs, namelen, myrank
    integer, dimension(MPI_STATUS_SIZE) :: stat
    integer :: tag = 1  ! tag for sending the data around
    integer :: nbeads   !number of beads per process
    integer :: procid   !ID to which we send or from which we recieve
    double precision, allocatable, dimension(:,:,:) :: buffer               ! buffer for tranfer of data between processes 

    CALL MPI_Comm_size ( MPI_COMM_WORLD, nprocs, error)
    CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)

    nbeads = n          !Number of beads this process has to take care of
    if (useMPI) then
        nbeads = n / nprocs
        allocate(buffer(d*p,d*p,nbeads))

        if (myrank .ne. root) then
                return  !should only be run by root
        end if
    endif

    ! =======================================================================
    ! SETTING UP PROGRESS INDICATORS
    ! =======================================================================
    
    totalentries = ((d*p)+1)*(d*p*n)
    entriestocalc = (d*p*n*((d*p)+1))/2
    WRITE(*,'(A,I12)') 'Number of entries in Hessian array: ', totalentries
    WRITE(99,'(A,I12)') 'Number of entries in Hessian array: ', totalentries
    WRITE(*,'(A,I12)') 'Of these, we now calculate: ', entriestocalc
    WRITE(99,'(A,I12)') 'Of these, we now calculate: ', entriestocalc

    ! Set up values to test when to write output to terminal and file:
    terminalreal = REAL(terminalint)
    fileint = 2*d*p*n
    filereal = REAL(fileint)

    entrycount = 0

    hessian = 0.0d0
    DO beadi = 1, nbeads
        CALL hessext(xmin(:,:,beadi), H)

        !debugging files:
        !------------------------------------------------------------------
        ! Write( hessfile, '(a,i10)' )  "Hessian",beadi
        ! OPEN (UNIT=2, FILE=hessfile, STATUS='REPLACE', ACTION='WRITE' )
        ! WRITE(2,*) H
        ! CLOSE (UNIT=2)
        !------------------------------------------------------------

        DO j = 1, p*d
            DO k = j, p*d
                hessian(p*d+2-j,p*d*(beadi-1)+k) = H(k-j+1,k)
                !diagonal is in last column
                !side diagonal is in last column but one
                ! ...

                ! Reporting progress to user:
                if (useMPI) then
                        ! +nprocs, because nprocs run in parallel 
                        ! => assume they are equally quick
                        entrycount = entrycount + nprocs
                else
                        entrycount = entrycount + 1
                end if
                percentcomplete = (REAL(entrycount)/REAL(entriestocalc))*100.0

                IF ( entrycount/terminalint == entrycount/terminalreal ) THEN
                    WRITE(*,2223) entrycount, entriestocalc-entrycount, percentcomplete
                    2223 FORMAT ('\rEntry of Hessian: ', I12, '; remaining: ', I12, '; (', F5.1, '% complete)', $ )
                    ! Dollar sign at end of previous command is included to suppress the automatic new line
                END IF

                IF ( entrycount/fileint == entrycount/filereal ) THEN
                    WRITE(99,2224) entrycount, entriestocalc-entrycount, percentcomplete
                    2224 FORMAT ('Entry of Hessian: ', I12, '; remaining: ', I12, '; (', F5.1, '% complete)' )
                END IF

            END DO
        END DO
    END DO

    !collect work from other processes:
    if (nprocs .gt. 1 .and. useMPI) then
         do procid=1,nprocs-1        !ie this loop is done nproc-1 times!
              call MPI_RECV(buffer, nbeads*d*p*d*p, MPI_DOUBLE_PRECISION,procid,tag,MPI_COMM_WORLD,stat,error)
              !copy to apropriate place in hessian matrix:
              do ib=1,nbeads
                   beadi=procid*nbeads+ib
                   DO j = 1, p*d
                       DO k = j, p*d
                           hessian(p*d+2-j,p*d*(beadi-1)+k) = buffer(k-j+1,k,ib)
                       enddo !k
                   enddo !j
                   
                !debugging files:
                !------------------------------------------------------------------
                ! DO j = 1, p*d
                !         DO k = 1, p*d
                !                 H(j,k) = buffer(j,k,ib)
                !         enddo !k
                ! enddo !j
                ! Write( hessfile, '(a,i10)' )  "Hessian",beadi
                ! OPEN (UNIT=2, FILE=hessfile, STATUS='REPLACE', ACTION='WRITE' )
                ! WRITE(2,*) H
                ! CLOSE (UNIT=2)
                !--------------------------------------------------------------

              enddo !ib
         enddo !procid
    endif

    WRITE(*,*)
    WRITE(*,'(A)') 'Remaining Hessian entries (-1s and +2s) will now be completed'
    WRITE(99,'(A)') 'Remaining Hessian entries (-1s and +2s) will now be completed'

    IF (n > 1) THEN
        hessian(1,:) = - 1.0d0 / (betan*hbar)**2
        hessian(d*p+1,:) = hessian(d*p+1,:) + 2.0d0 / (betan*hbar)**2
    END IF

    WRITE(*,'(A)') 'HESSIAN COMPLETE'
    WRITE(99,'(A)') 'HESSIAN COMPLETE'

    END SUBROUTINE hesscalc

! ===========================================================================================================

    SUBROUTINE hesscalc_child ( xmin )
    USE funcs						! Miscellaneous functions module
    USE grads						! Gradient calculation module
    USE mpi
    use parameters, only: root
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! Calculate the hessian matrix for a selected number of beads and send them
    ! to the main subroutine above
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  24 Nov 2012    M. Herbst      Original subroutine
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Calling variables:
    double precision,DIMENSION(d,p,n),INTENT(IN) :: xmin			! Minimum geometry
    
    ! Local variables (user may wish to alter the assigned values):
    double precision :: delta = 0.0001					! Finite difference in co-ordinate

    ! Local variables (no assigned values):
    INTEGER :: ib,j, k							! Indices to loop over to create banded Hessian
    double precision,DIMENSION(d*p,d*p) :: H				! Hessian of one bead

    !MPI setup
    INTEGER :: error, nprocs, namelen, myrank
    integer, dimension(MPI_STATUS_SIZE) :: stat
    integer :: tag = 1  ! tag for sending the data around
    integer :: nbeads   !number of beads per process
    integer :: offset    !offset for this process
    double precision, allocatable, dimension(:,:,:) :: buffer               ! buffer for tranfer of data between processes 

    CALL MPI_Comm_size ( MPI_COMM_WORLD, nprocs, error)
    CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)

    if (myrank .eq. root) then
        ! the root process should not call this
        return
    endif

    if (root .ne. 0) then
            write (*,*) "The routine hesscalc_child assumes the id of the root process is zero which is not the case."
            stop
    end if

    nbeads = n / nprocs !No of beads this process has to take care of
    allocate(buffer(d*p,d*p,nbeads))
    offset = myrank*nbeads      !Note: myrank is 0based

    buffer = 0.0d0
    do ib = 1,nbeads
        CALL hessext(xmin(:,:,ib+offset), H)
        DO j = 1, p*d
            DO k = 1, p*d
                buffer(j,k,ib) = H(j,k)
            END DO
        END DO
    END DO

    !send work to root
    call MPI_SEND(buffer, nbeads*d*p*d*p,MPI_DOUBLE_PRECISION,root,tag,MPI_COMM_WORLD,error)
    END SUBROUTINE hesscalc_child

! ===========================================================================================================

    SUBROUTINE hesscalcsingle ( xsinglemin, hessian )
    USE funcs						! Miscellaneous functions module
    USE grads						! Gradient calculation module
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To take an an xyz file representing a linear polymer and calculate the
    ! Hessian matrix of second derivatives.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  7 Nov 2011     A. Reid        Development of subroutine
    !  15 Nov 2011    A. Reid        Completion of initial version
    !  20 Jan 2012    A. Reid        Updated: shared variables / non-scaled x
    !  31 Jan 2012    A. Reid        Added debugging code
    !  27 Feb 2012    A. Reid        Grad calculations moved to subroutines
    !  31 Mar 2012    J. Richardson  Simplified
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Calling variables:
    double precision,DIMENSION(d*p),INTENT(IN) :: xsinglemin		! Minimum geometry
    double precision,DIMENSION((d*p),(d*p)),INTENT(INOUT) :: hessian	! fixed-ended mass-weighted linear polymer hessian (upper) (band storage)
    
    ! Local variables (user may wish to alter the assigned values):
    !CHARACTER(len=10) :: terminalwrite = '1'				! Interval between writes to terminal
    CHARACTER(len=10) :: filewrite = '1000'				! Interval between writes to file
    double precision :: delta = 1.0d-3					! Finite difference in co-ordinate

    ! Local variables (no assigned values):
    INTEGER :: beadi							! Index of time slices
    INTEGER :: j, k							! Indices to loop over to create banded Hessian
    INTEGER :: totalentries						! Total no. of entries in Hessian
    double precision :: percentcomplete					! Percentage completion
    double precision :: totalentriesreal				! Total no. of entries in Hessian (as a real number)
    double precision :: entrycountreal					! Hessian entry no. (as a real number)
    double precision,DIMENSION(d*p,d*p) :: H				! Hessian of one bead


    ! =======================================================================
    ! SETTING UP PROGRESS INDICATORS
    ! =======================================================================
    
    totalentries = (d*p)*(d*p)
    WRITE(*,'(A,I5,A)') 'Hessian array has ', totalentries, ' entries'
    WRITE(*,'(A)') 'Calculating: please wait...'
    WRITE(99,'(A,I5,A)') 'Hessian array has ', totalentries, ' entries'
    WRITE(99,'(A)') '-----------------------------------------------------------------------'

    ! Set up values to test when to write output to terminal and file:
    
    CALL hessext(xsinglemin, hessian)

    END SUBROUTINE hesscalcsingle

! ===========================================================================================================

    SUBROUTINE hessext ( x, H )
    !USE pes_shell
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To calculate the second derivative of the PES at one bead
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  24 Feb 2012    A. Reid        Development of subroutine
    !  31 Mar 2012    J. Richardson  Modified to apply to just one bead
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Calling variables:
    double precision,DIMENSION(d*p),INTENT(IN) :: x        ! Geometry of a bead
    double precision,DIMENSION(d*p) :: xpes        ! Geometry of a bead in pes coords
    double precision,DIMENSION(d*p,d*p),INTENT(OUT) :: H   ! Hessian of one bead
    
    ! Local variables:
    INTEGER :: dimi
    INTEGER :: atomi
    INTEGER :: dimj
    INTEGER :: atomj

    xpes = x
    call pes_obj%vec_to_pes(xpes)	! transforms vector from coordinate convention used by program to coordinate convention used by the pes potential
    call pes_obj%hess(xpes, hesseps, H)	! hesseps comes from command line, or parameters.f90 if command line arg not given
    call pes_obj%mat_from_pes(H)	! transforms Hessian from coord convention used by pes to coord convention used by program
    ! mass-weight
    FORALL(atomi=1:p,dimi=1:d,atomj=1:p,dimj=1:d) &
            H((atomi-1)*d+dimi,(atomj-1)*d+dimj) = H((atomi-1)*d+dimi,(atomj-1)*d+dimj) / sqrt(masses(atomi)*masses(atomj))

    END SUBROUTINE hessext

! ===========================================================================================================

    SUBROUTINE hessout ( d, p, n, hessian, hessfile )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To take a 3D array representing a linear-polymer and output an xyz file
    ! that can be read by VMD.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  4 Nov 2011     A. Reid        Development of subroutine
    !  4 Nov 2011     A. Reid        Amended to take array of different shape
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Calling variables:
    INTEGER,INTENT(IN) :: d					! No. of dimensions
    INTEGER,INTENT(IN) :: n					! No. of time slices / f-dimensional beads
    INTEGER,INTENT(IN) :: p					! No. of atoms
    CHARACTER(len=32),INTENT(IN) :: hessfile			! File to create
    !double precision,DIMENSION(d*p*n,d*p*n),INTENT(IN) :: hessian	! Hessian of final linear-polymer
    double precision,DIMENSION(((d*p)+1),d*p*n),INTENT(IN) :: hessian	! Hessian of final linear-polymer
    
    ! Local variables:
    CHARACTER(80) :: msg					! Status message
    INTEGER :: status						! IO status message
    
    
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    WRITE(99,'(A)') 'Writing the hessian to the output file'
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    
    OPEN (UNIT=13, FILE=hessfile, STATUS='REPLACE', ACTION='WRITE', IOSTAT=status, IOMSG=msg )
    
        WRITE(13,*) hessian
    
    CLOSE (UNIT=13)
    
    END SUBROUTINE hessout

! ===========================================================================================================

    SUBROUTINE diag ( d, p, no, matrix, evectors, evalues )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! Uses the ACML library to diagonalize a matrix.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  ?? ??? ????    M. Cvitas      Original subroutine developed
    !  15 Nov 2011    A. Reid        Adaptation of Marko's subroutine
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Calling variables:
    INTEGER,INTENT(IN) :: d
    INTEGER,INTENT(IN) :: p
    INTEGER,INTENT(IN) :: no
    double precision,DIMENSION(d*p*no, d*p*no),INTENT(IN) :: matrix	! Matrix to be diagonalized
    double precision,DIMENSION(d*p*no, d*p*no),INTENT(OUT) :: evectors	! Computed eigenvectors of matrix
    double precision,DIMENSION(SIZE(matrix,1)),INTENT(OUT) :: evalues	! Computed eigenvalues of matrix
    
    ! Local variables:
    INTEGER :: n
    CHARACTER(1) :: jobz, uplo
    double precision,DIMENSION(1+6*SIZE(matrix,1)+2*SIZE(matrix,1)**2) :: work
    INTEGER,DIMENSION(3+5*SIZE(matrix,1)) :: iwork
    INTEGER :: lwork, liwork
    INTEGER :: lda, info
    
    ! copy the diagonal and off-diagonal elements
    jobz='V'
    uplo='U'
    
    n = SIZE(matrix,1)
    lda = n
    
    lwork = 1+6*n+2*n*n
    liwork = 3+5*n
    
    ! copy the matrix
    evectors = matrix
    
    ! lapack call
    CALL DSYEVD(jobz,uplo,n,evectors,lda,evalues,work,lwork,iwork,liwork,info)
    
    ! results check
    IF (info.NE.0) THEN
        WRITE(*,'(A,I10)') 'symmetric matrix diagonalization error: info =', info
        WRITE(99,'(A,I10)') 'symmetric matrix diagonalization error: info =', info
        STOP
    END IF
    
    END SUBROUTINE diag

! ===========================================================================================================

    SUBROUTINE diagband ( ab, w, z, ldz )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! Uses the ACML library to diagonalize a banded matrix, denoted A and
    ! stored in an array in banded form.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  28 Feb 2012    A. Reid        Development of subroutine from 'diag'
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================

    ! Local variables (1/3):
    CHARACTER(1) :: jobz = 'N'					! (1) 'N' = compute evalues only; 'V' = compute evalues and evectors
    CHARACTER(1) :: uplo = 'U'					! (2) 'U' = upper triangle of A is stored; 'L' = lower triangle of A is stored
    INTEGER :: order						! (3) order = the order of the matrix A ('n' not used, to avoid confusion with beads)
    INTEGER :: kd						! (4) kd = number of non-zero super-diagonals (if uplo = 'U') or sub-diagonals (if uplo = 'L')
    
    ! Calling variables (1/2):
    double precision,DIMENSION((d*p)+1,d*p*n),INTENT(IN) :: ab	! (5) ab = banded array representing the matrix to be diagonalized

    ! Local variables (2/3):
    INTEGER :: ldab						! (6) ldab = leading dimension of the array AB

    ! Calling variables (2/2):
    double precision,DIMENSION(d*p*n),INTENT(OUT) :: w		! (7) Computed eigenvalues of matrix, in ascending order
    double precision,DIMENSION(ldz,d*p*n),INTENT(OUT) :: z	! (8) Computed eigenvectors of matrix (variable size, depending on ldz)
    INTEGER,INTENT(IN) :: ldz					! (9) Number of eigenvectors to output in array 'z'
    
    ! Local variables (3/3):
    double precision,DIMENSION( (3*(d*p*n)) - 2 ) :: work	! (10) work = workspace array for lapack usage
    INTEGER :: info						! (11) when info = 0, diagonalization has exited successfully


    order = d*p*n
    kd = d*p               
    ldab = kd+1							! ldab = leading dimension of array AB; hence, ldab = kd+1
    
    ! lapack call
    CALL DSBEV(jobz,uplo,order,kd,ab,ldab,w,z,ldz,work,info)
    
    ! results check
    IF (info.NE.0) THEN
        WRITE(*,'(A,I10)') 'symmetric matrix diagonalization error: info =', info
        WRITE(99,'(A,I10)') 'symmetric matrix diagonalization error: info =', info
        STOP
    END IF
    
    END SUBROUTINE diagband

! ===========================================================================================================

! IN DEVELOPMENT

    SUBROUTINE diagsymm ( a, w, z )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! Uses the ACML library to diagonalize a banded matrix, denoted A and
    ! stored in an array in banded form.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  28 Feb 2012    A. Reid        Development of subroutine from 'diag'
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================

    ! Local variables (1/3):
    CHARACTER(1) :: jobz = 'N'					! (1) 'N' = compute evalues only; 'V' = compute evalues and evectors
    CHARACTER(1) :: uplo = 'U'					! (2) 'U' = upper triangle of A is stored; 'L' = lower triangle of A is stored
    INTEGER :: order						! (3) order = the order of the matrix A ('n' not used, to avoid confusion with beads)
    
    ! Calling variables (1/3):
    double precision,DIMENSION((d*p),(d*p)),INTENT(IN) :: a	! (4) a = banded array representing the matrix to be diagonalized. On exit, if JOBZ = 'V',
								! then if INFO = 0, A contains the orthonormal eigenvectors of the matrix A.
    ! Local variables (2/3):
    INTEGER :: lda						! (5) lda = leading dimension of the array A

    ! Calling variables (2/3):
    double precision,DIMENSION(d*p),INTENT(OUT) :: w		! (6) Computed eigenvalues of matrix, in ascending order
    
    ! Local variables (3/3):
    double precision,DIMENSION( (3*(d*p)) ) :: work		! (7) work = workspace array for lapack usage
    INTEGER :: lwork						! (8) lwork = length of the work array; lwork >= MAX(1,3*order-1)
    INTEGER :: info						! (9) when info = 0, diagonalization has exited successfully

    ! Calling variables (3/3):
    double precision,DIMENSION((d*p),(d*p)),INTENT(OUT) :: z	! Copy of a, if eigenvectors are calculated

    order = d*p
    lda = order							! lda = leading dimension of array A; hence, lda = order
    lwork = 3*(d*p)

    ! lapack call
    CALL DSYEV(jobz,uplo,order,a,lda,w,work,lwork,info)
    
    ! results check
    IF (info.NE.0) THEN
        WRITE(*,'(A,I10)') 'symmetric matrix diagonalization error: info =', info
        WRITE(99,'(A,I10)') 'symmetric matrix diagonalization error: info =', info
        STOP
    END IF

    IF (jobz == 'V') THEN
        z = a							! Copy a, which now contains eigenvectors
    END IF
    
    END SUBROUTINE diagsymm

! ===========================================================================================================

END MODULE hess
