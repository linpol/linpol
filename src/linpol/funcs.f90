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
MODULE funcs
USE parameters
USE prep
IMPLICIT NONE
private
public :: calcpos,onetothree,onetotwo,threetoonequick,threetoone,twotoone,&
        xscaling,xyzout,startendxyzout,getdatetime,outputgeometry, dateToString
!
! =======================================================================
! PURPOSE OF THIS MODULE
! =======================================================================
! To contain the functions and subroutines to accomplish miscellaneous
! tasks.
! 
! Contents:
!     calcpos
!     onetothree
!     onetotwo
!     threetoone
!     threetoonequick
!     xscaling
!     xyzout
!     calcdatetime
!     outputgeometry
!
! Record of revisions:
!
!  Date           Person         Change implemented
!  ----           ------         ------------------
!  21 Nov 2011    A. Reid        Development of module
!  29 Nov 2011    A. Reid        Conversion to use variables from modules
!  20 Dec 2011    A. Reid        Addition of calcdatetime subroutine
!
! =======================================================================

! =======================================================================
! DECLARATION OF GLOBAL VARIABLES AVAILABLE TO OTHER PROGRAM UNITS
! =======================================================================

 CHARACTER(len=15),save :: datetime = ''

CONTAINS

! ===========================================================================================================

    FUNCTION calcpos ( d, p, dim, atom, bead )
    !
    ! =======================================================================
    ! PURPOSE OF THIS FUNCTION
    ! =======================================================================
    ! To calculate the position in a 1D array where the value of a particular
    ! degree of freedom is stored.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  11 Nov 2011    A. Reid        Development of function
    !
    ! =======================================================================

    ! Calling parameters:
    INTEGER,INTENT(IN) :: d                                ! No. of dimensions
    INTEGER,INTENT(IN) :: p                                ! No. of atoms
    INTEGER,INTENT(IN) :: atom                             ! Index of atoms
    INTEGER,INTENT(IN) :: dim                              ! Index of dimensions
    INTEGER,INTENT(IN) :: bead                             ! Index of time slices
    INTEGER :: calcpos                                     ! Returned position

    calcpos = (d*(atom-1) + dim) + ((bead-1)*d*p)
    ! First term: calculates how far along the bead the position is
    ! Second term: calculates how many positions to add on if not in the first bead

    END FUNCTION calcpos

! ===========================================================================================================

    SUBROUTINE onetothree ( d, p, n, xone, xthree )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To take in a 1D array with progressive bead, atom and dimension
    ! subunits (see example below) and convert it to a 3D array indexed by
    ! dimension, atom and bead (also below).
    !
    ! Example: 2 dimensions, 2 atoms, 2 beads
    ! 1D: d1a1b1,d2a1b1,d1a2b1,d2a2b1,d1a1b2,d2a1b2,d1a2b2,d2a2b2
    ! 3D: a 2x2x2 array that can be visualized as a cube containing 8 entries
    !
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  18 Nov 2011    A. Reid        Development of subroutine
    !
    ! =======================================================================

    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================

    ! Calling variables:
    INTEGER,INTENT(IN) :: d                               ! No. of dimensions
    INTEGER,INTENT(IN) :: p                               ! No. of atoms
    INTEGER,INTENT(IN) :: n                               ! No. of time slices / f-dimensional beads
    double precision,DIMENSION(d*p*n),INTENT(IN) :: xone              ! Ring polymer geometry in a 1D array
    double precision,DIMENSION(d,p,n),INTENT(OUT) :: xthree           ! Ring polymer geometry in a 3D array

    ! Local variables:
    INTEGER :: atom                                       ! Index of atoms
    INTEGER :: dim                                        ! Index of dimensions
    INTEGER :: bead                                       ! Index of time slices

    beadloop: DO bead = 1, n
        atomloop: DO atom = 1, p
            dimloop: DO dim = 1, d

                xthree(dim,atom,bead) = xone( calcpos(d,p,dim,atom,bead) )

            END DO dimloop
        END DO atomloop
    END DO beadloop

    END SUBROUTINE onetothree

! ===========================================================================================================

    SUBROUTINE onetotwo ( d, p, xone, xtwo )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To take in a 1D array with progressive atom and dimension subunits (see
    ! example below) and convert it to a 2D array indexed by dimension and
    ! atom(also below).
    !
    ! Example: 2 dimensions, 2 atoms, 2 beads
    ! 1D: d1a1,d2a1,d1a2,d2a2
    ! 2D: a 2x2 array that can be visualized as a square containing 4 entries
    !
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  5 Dec 2011     A. Reid        Development of subroutine
    !
    ! =======================================================================

    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================

    ! Calling variables:
    INTEGER,INTENT(IN) :: d                               ! No. of dimensions
    INTEGER,INTENT(IN) :: p                               ! No. of atoms
    double precision,DIMENSION(d*p),INTENT(IN) :: xone              ! Ring polymer geometry in a 1D array
    double precision,DIMENSION(d,p),INTENT(OUT) :: xtwo           ! Ring polymer geometry in a 3D array

    ! Local variables:
    INTEGER :: atom                                       ! Index of atoms
    INTEGER :: dim                                        ! Index of dimensions

    atomloop: DO atom = 1, p
        dimloop: DO dim = 1, d

            xtwo(dim,atom) = xone( calcpos(d,p,dim,atom,1) )

        END DO dimloop
    END DO atomloop

    END SUBROUTINE onetotwo

! ===========================================================================================================

    SUBROUTINE threetoone ( d, p, n, xthree, xone )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To take in a 3D array indexed by dimension, atom and bead (also below)
    ! and convert it to 1D array with progressive bead, atom and dimension
    ! subunits (see example below)
    ! 
    ! Example: 2 dimensions, 2 atoms, 2 beads
    ! 3D: a 2x2x2 array that can be visualized as a cube containing 8 entries
    ! 1D: d1a1b1,d2a1b1,d1a2b1,d2a2b1,d1a1b2,d2a1b2,d1a2b2,d2a2b2
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
    INTEGER,INTENT(IN) :: d                               ! No. of dimensions
    INTEGER,INTENT(IN) :: p                               ! No. of atoms
    INTEGER,INTENT(IN) :: n                               ! No. of time slices / f-dimensional beads
    double precision,DIMENSION(d,p,n),INTENT(IN) :: xthree            ! Ring polymer geometry in a 3D array
    double precision,DIMENSION(d*p*n),INTENT(OUT) :: xone             ! Ring polymer geometry in a 1D array

    ! Local variables:
    INTEGER :: atom                                       ! Index of atoms
    INTEGER :: dim                                        ! Index of dimensions
    INTEGER :: bead                                       ! Index of time slices

    beadloop: DO bead = 1, n
        atomloop: DO atom = 1, p
            dimloop: DO dim = 1, d

                xone( calcpos(d,p,dim,atom,bead) ) = xthree(dim,atom,bead)

            END DO dimloop
        END DO atomloop
    END DO beadloop

    END SUBROUTINE threetoone

! ===========================================================================================================

    SUBROUTINE twotoone ( d, p, xtwo, xone )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To take in a 2D array indexed by dimension and atom (also below)
    ! and convert it to 1D array with progressive atom and dimension
    ! subunits (see example below)
    ! 
    ! Example: 2 dimensions, 2 atoms
    ! 2D: a 2x2 array that can be visualized as a square containing 4 entries
    ! 1D: d1a1,d2a1,d1a2,d2a2
    !
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  5 Dec 2011     A. Reid        Development of subroutine
    !
    ! =======================================================================

    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================

    ! Calling variables:
    INTEGER,INTENT(IN) :: d                               ! No. of dimensions
    INTEGER,INTENT(IN) :: p                               ! No. of atoms
    double precision,DIMENSION(d,p),INTENT(IN) :: xtwo               ! Cluster geometry in a 2D array
    double precision,DIMENSION(d*p),INTENT(OUT) :: xone             ! Cluster geometry in a 1D array

    ! Local variables:
    INTEGER :: atom                                       ! Index of atoms
    INTEGER :: dim                                        ! Index of dimensions

    atomloop: DO atom = 1, p
        dimloop: DO dim = 1, d

            xone( calcpos(d,p,dim,atom,1) ) = xtwo(dim,atom)

        END DO dimloop
    END DO atomloop

    END SUBROUTINE twotoone

! ===========================================================================================================

    SUBROUTINE threetoonequick ( d, p, n, xthree, xone )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To take in a 3D array indexed by dimension, atom and bead (also below)
    ! and convert it to 1D array with progressive bead, atom and dimension
    ! subunits (see example below)
    ! 
    ! Example: 2 dimensions, 2 atoms, 2 beads
    ! 3D: a 2x2x2 array that can be visualized as a cube containing 8 entries
    ! 1D: d1a1b1,d2a1b1,d1a2b1,d2a2b1,d1a1b2,d2a1b2,d1a2b2,d2a2b2
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
    INTEGER,INTENT(IN) :: d                               ! No. of dimensions
    INTEGER,INTENT(IN) :: p                               ! No. of atoms
    INTEGER,INTENT(IN) :: n                               ! No. of time slices / f-dimensional beads
    double precision,DIMENSION(d*p*n),INTENT(IN) :: xthree            ! Ring polymer geometry in a 3D array
    double precision,DIMENSION(d*p*n),INTENT(OUT) :: xone             ! Ring polymer geometry in a 1D array

    xone = xthree

    END SUBROUTINE threetoonequick

! ===========================================================================================================

!    SUBROUTINE xscaling ( x, n, d, p, masses, xscaled )
    SUBROUTINE xscaling ( x, xscaled )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To scale the values of the 3D array x in order to calculate the
    ! 'internal' part of the linear-polymer potential.
    ! 
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  1 Nov 2011     A. Reid        Development of subroutine
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Calling variables:
!    INTEGER,INTENT(IN) :: n                               ! No. of time slices / f-dimensional beads
!    INTEGER,INTENT(IN) :: d                               ! No. of dimensions
!    INTEGER,INTENT(IN) :: p                               ! No. of atoms
double precision,DIMENSION(d,p,n),INTENT(IN) :: x                 ! Unscaled linear-polymer geometry
!    REAL(kind=dp),DIMENSION(p),INTENT(IN) :: masses                ! Type of element
double precision,DIMENSION(d,p,n),INTENT(OUT) :: xscaled          ! Mass- and unit-scaled linear-polymer geometry
    
    ! Local variables:
    INTEGER :: ind                                        ! Index of DO loops
    INTEGER :: indsub                                     ! Index of sub DO loops
    INTEGER :: indsubsub                                  ! Index of sub-sub DO loops
    double precision :: massroot                                      ! Square root of the mass of the degree of freedom being dealt with
    
    
    
    xscaled = x
    
    massscaling: DO ind = 1, p                             ! DO loop from 1 to p (i.e. over atoms), each of which has a different mass
        massroot = SQRT( masses( ind ) )
    
        massscalingsub: DO indsub = 1, d                      ! DO loop from 1 to d (i.e. over dimensions)
            massscalingsubsub: DO indsubsub = 1, n                  ! DO loop from 1 to n (i.e. over beads)
                xscaled( indsub, ind, indsubsub ) =  massroot * xscaled( indsub, ind, indsubsub )     ! Scale by sqrt of mass
            END DO massscalingsubsub
        END DO massscalingsub
    END DO massscaling
    
    END SUBROUTINE xscaling
    
! ===========================================================================================================

    SUBROUTINE xyzout ( x , n , d , p, beta, atomtags, outputfile )
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
    !
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Calling variables:
    INTEGER,INTENT(IN) :: d                               ! No. of dimensions
    INTEGER,INTENT(IN) :: n                              ! No. of time slices / f-dimensional beads
    INTEGER,INTENT(IN) :: p                              ! No. of atoms
    double precision,DIMENSION(d,p,n),INTENT(IN) :: x             ! Linear-polymer with intermediate beads added
    double precision,dimension(d,p) :: bead
    CHARACTER(len=32),INTENT(IN) :: outputfile            ! File to create
    double precision,INTENT(IN) :: beta                              ! Beta
    CHARACTER(len=2),DIMENSION(p),INTENT(IN) :: atomtags  ! Atom tags
    
    ! Local variables:
    INTEGER :: atom                                          ! Index of atoms
    INTEGER :: dim                                          ! Index of dimensions
    INTEGER :: beadi                                         ! Index of time slices
    CHARACTER(80) :: msg                                  ! Status message
    INTEGER :: status                                     ! IO status message
    
    
    OPEN (UNIT=2, FILE=outputfile, STATUS='REPLACE', ACTION='WRITE', IOSTAT=status, IOMSG=msg )
    
    !if (transformed_to_pot) then
    !        call pes_obj%vec_from_pes(atomtags) !transform forward for writing
    !end if
    
    beadloop: DO beadi = 1,n
    
        WRITE(2,'(I0)') p
        WRITE(2,'(F0.3)') beta

        bead = x(:,:,beadi)
        if (transformed_to_pot) then
                call pes_obj%vec_from_pes( bead )
        end if

        atomloop: DO atom = 1,p
            WRITE(2,*) atomtags(atom), bead(:,atom)*auang
        END DO atomloop
    END DO beadloop
    
    CLOSE (UNIT=2)
    
    !if (transformed_to_pot) then
    !        call pes_obj%vec_to_pes(atomtags)   !transform back to what it was
    !end if
    END SUBROUTINE xyzout

! ===========================================================================================================

    SUBROUTINE startendxyzout ( x, d, p, beta, atomtags, outputfile )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To take a 2D array representing a starting or ending geometry and
    ! output an xyz file that can be read by VMD.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  5 Dec 2011     A. Reid        Development of subroutine
    !
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================
    
    ! Calling variables:
    INTEGER,INTENT(IN) :: d                               ! No. of dimensions
    INTEGER,INTENT(IN) :: p                              ! No. of atoms
    double precision,DIMENSION(d,p),INTENT(IN) :: x        ! Starting or ending geometry
    CHARACTER(len=32),INTENT(IN) :: outputfile            ! File to create
    CHARACTER(len=2),DIMENSION(p),INTENT(IN) :: atomtags  ! Atom tags
    double precision,INTENT(IN) :: beta                              ! Beta

    call xyzout (x,1,d,p,beta,atomtags,outputfile)


    
!     ! Local variables:
!     INTEGER :: atom                                          ! Index of atoms
!     INTEGER :: dim                                          ! Index of dimensions
!     CHARACTER(80) :: msg                                  ! Status message
!     INTEGER :: status                                     ! IO status message
!     
!     OPEN (UNIT=3, FILE=outputfile, STATUS='REPLACE', ACTION='WRITE', IOSTAT=status, IOMSG=msg )
!     WRITE(3,'(I0)') p
!     WRITE(3,'(F0.3)') beta
! 
!     if (transformed_to_pot) then
!         call pes_obj%vec_from_pes( x )
!     end if
!     atomloop: DO atom = 1,p
!         WRITE(3,*) atomtags(atom), x(:,atom)*auang
!     END DO atomloop
!     
!     CLOSE (UNIT=3)
    
    END SUBROUTINE startendxyzout

! ===========================================================================================================

        function getdatetime()
                !
                ! =======================================================================
                ! PURPOSE OF THIS FUNCTION
                ! =======================================================================
                ! Get a unique date and time. Make sure that each time function is
                ! called the same date/time is returned (the one of the initial call)
                ! Intention is the name of the output files.
                !
                ! Record of revisions:
                !
                !  Date           Person         Change implemented
                !  ----           ------         ------------------
                !  23 Dec 2011    M. Herbst      Development of function
                !
                ! =======================================================================
                CHARACTER(len=8)  :: date
                CHARACTER(len=10) :: time
                CHARACTER(len=15) :: getdatetime

                if (datetime == '') then
                        CALL date_and_time(DATE=date,TIME=time)
                        datetime = date//'-'//time(1:6)
                end if
                getdatetime = datetime
        end function

! ===========================================================================================================

        function dateToString(seconds) result(str)
                double precision ,intent(in) :: seconds         ! the number of seconds representing the date we want to write.
                INTEGER :: days, hours, mins, secs  ! For expressing the program duration in integers
                character(len=46) :: Str

                days = INT(seconds/(60d0*60d0*24d0))
                hours = INT(seconds/(60d0*60d0)) - (days*24d0)
                mins = INT(seconds/(60d0)) - (days*24d0*60d0) - (hours*60d0)
                secs = INT(seconds) - (days*24d0*60d0*60d0) - (hours*60d0*60d0) - (mins*60d0)

                if ((mins == 0 .and. hours == 0) .and. (days == 0)) then
                        WRITE(Str,'(35x,i2,a)') secs, ' sec(s)'
                else if ((hours == 0) .and. (days == 0)) then
                        write(Str,'(24x,I2,a,i2,a)') mins, ' min(s), ', secs, ' sec(s)'
                else if (days == 0) then
                        WRITE(Str,'(12x,I2,a,I2,a,i2,a)') hours, ' hour(s), ', mins, &
                                ' min(s), ', secs, ' sec(s)'
                else
                        WRITE(Str,'(I3,A,I2,A,I2,A,I2,A)') days, ' day(s), ', hours, &
                                ' hour(s), ', mins, ' min(s), ', secs, ' sec(s)'
                end if

                !IF ((days == 0) .AND. (hours == 0) .AND. (mins == 0)) THEN
                !    WRITE(Str,'(I2,A)') secs, ' sec(s)'
                !ELSEIF ((days == 0) .AND. (hours == 0)) THEN
                !    WRITE(Str,'(I2,A,I2,A)') mins, ' min(s) ', secs, ' sec(s)'
                !ELSEIF (days == 0) THEN
                !    WRITE(Str,'(I2,A,I2,A,I2,A)') hours, ' hour(s), ', mins, ' min(s) ', secs, ' second(s)'
                !ELSE
                !    WRITE(Str,'(I3,A,I2,A,I2,A,I2,A)') days, ' day(s), ', hours, ' hour(s), ', mins, ' min(s) ', secs, ' sec(s)'
                !END IF
        end function

! ===========================================================================================================

    SUBROUTINE outputgeometry ( xmin , n , d , p, outputgeometryfile )
    !
    ! =======================================================================
    ! PURPOSE OF THIS SUBROUTINE
    ! =======================================================================
    ! To output the optimized geometry to a text file that can be used as the
    ! starting point for future optimizations with more beads.
    !
    ! Record of revisions:
    !
    !  Date           Person         Change implemented
    !  ----           ------         ------------------
    !  21 Dec 2011    A. Reid        Development of subroutine
    !
    ! =======================================================================
    
    ! =======================================================================
    ! DECLARATION OF VARIABLES
    ! =======================================================================

    ! Calling variables:
    INTEGER,INTENT(IN) :: d                              ! No. of dimensions
    INTEGER,INTENT(IN) :: n                              ! No. of time slices / f-dimensional beads
    INTEGER,INTENT(IN) :: p                              ! No. of atoms
    double precision,DIMENSION(d,p,n),INTENT(IN) :: xmin    ! Linear-polymer as optimized by L-BFGS

    CHARACTER(len=64),INTENT(IN) :: outputgeometryfile

    ! Local variables:
    INTEGER :: atom                                       ! Index of atoms
    INTEGER :: dim                                        ! Index of dimensions
    INTEGER :: bead                                       ! Index of time slices
    CHARACTER(80) :: msg                                  ! Status message
    INTEGER :: status                                     ! IO status message
    
    
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    WRITE(99,'(A)') 'Creating geometry output file for use in future optimizations'
    WRITE(99,'(A)') '-----------------------------------------------------------------------'
    
    
    OPEN (UNIT=4, FILE=outputgeometryfile, STATUS='REPLACE', ACTION='WRITE', IOSTAT=status, IOMSG=msg )
    
    
    beadloop: DO bead = 1,n
        atomloop: DO atom = 1,p

            WRITE(4,*) xmin(1,atom,bead), xmin(2,atom,bead), xmin(3,atom,bead)

        END DO atomloop
    END DO beadloop
    
    CLOSE (UNIT=4)


    END SUBROUTINE outputgeometry

! ===========================================================================================================

END MODULE funcs
