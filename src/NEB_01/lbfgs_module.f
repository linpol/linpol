      module lbfgs_module
      use potential_optim_module, only: pot_and_grad, set_min_pot
      use opt_events ! module to connect to this one and the main program

      implicit none

      private
      public :: optimize, optim_init

*** LBFGS PARAMETERS
      integer, save :: msave, ndim
      double precision, save :: eps
      integer, save :: MAXITER

      double precision, save :: DGUESS=1.0D0
      double precision, save :: MAXMOVE=10.0D0
      double precision, save :: MAXERISE=1.0D0

      contains

****************************************************************************
*** THIS IS THE DRIVER SUBROUTINE FOR LBFGS by Nocedal
****************************************************************************
      subroutine optim_init(ichoice,m,epsilon,maxiterations)
      integer :: ichoice
      integer, intent(in), optional :: m
      double precision, intent(in), optional :: epsilon
      integer, intent(in), optional :: maxiterations

      call set_min_pot(ichoice)

      if (present(m)) then
         msave=m   
      else  
         msave=5
      end if

      if (present(epsilon)) then
         eps=epsilon
      else
         eps=1.0D-5
      end if

      if (present(maxiterations)) then
         MAXITER=maxiterations
      else
         MAXITER=10000
      end if
      end subroutine optim_init

****************************************************************************
      subroutine optimize(x,f,g,m,epsilon,maxiterations)
      use opt_events

      double precision, dimension(:) :: x, g
      double precision :: f

      integer, optional :: m
      double precision, optional :: epsilon
      integer, intent(in), optional :: maxiterations

      integer :: ifail

*** set lbfgs parameters optionally
      ndim=size(x)

*** transfer n = no of dimensions and m (saved steps: between 3 and 7)
      if (present(m)) then
         msave=m
      else
         msave=3
      end if

*** stopping criterion
      if (present(epsilon)) then
         eps=epsilon
      else
         eps=1.0D-5
      end if

      if (present(maxiterations)) then
         MAXITER=maxiterations
      else
         MAXITER=10000
      end if

      ifail=0
      call LBFGS(ndim,msave,x,f,g,eps,ifail)

      if (ifail.eq.1) then
              call opt_event_opt_restarted(ifail)
              call LBFGS(ndim,msave,x,f,g,eps,ifail)

              if (ifail.eq.2) then
                      call opt_event_opt_failed(ifail)
                      stop 1
              end if
      end if
      end subroutine optimize

****************************************************************************
C     ----------------------------------------------------------------------
C     This file contains the LBFGS algorithm and supporting routines
C
C     ****************
C     LBFGS SUBROUTINE
C     ****************
C
      SUBROUTINE LBFGS(N,M,X,F,G,EPS,IFAIL)

      INTEGER N,M,IFAIL
      DOUBLE PRECISION, INTENT(INOUT) :: X(N),G(N)
      DOUBLE PRECISION, ALLOCATABLE :: W(:), DIAG(:) ! W(N*(2*M+1)+2*M), DIAG(N)
      DOUBLE PRECISION F,EPS

C
C        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
C                          JORGE NOCEDAL
C                        *** July 1990 ***
C
C 
C     N       is an INTEGER variable that must be set by the user to the
C             number of variables. It is not altered by the routine.
C             Restriction: N>0.
C 
C     M       is an INTEGER variable that must be set by the user to
C             the number of corrections used in the BFGS update. It
C             is not altered by the routine. Values of M less than 3 are
C             not recommended; large values of M will result in excessive
C             computing time. 3<= M <=7 is recommended. Restriction: M>0.
C 
C     X       is a DOUBLE PRECISION array of length N. On initial entry
C             it must be set by the user to the values of the initial
C             estimate of the solution vector. On exit with IFLAG=0, it
C             contains the values of the variables at the best point
C             found (usually a solution).
C 
C     F       is a DOUBLE PRECISION variable. Before initial entry and on
C             a re-entry with IFLAG=1, it must be set by the user to
C             contain the value of the function F at the point X.
C 
C     G       is a DOUBLE PRECISION array of length N. Before initial
C             entry and on a re-entry with IFLAG=1, it must be set by
C             the user to contain the components of the gradient G at
C             the point X.
C 
C     EPS     is a positive DOUBLE PRECISION variable that must be set by
C             the user, and determines the accuracy with which the solution
C             is to be found. The subroutine terminates when
C
C                         ||G|| < EPS max(1,||X||),
C
C             where ||.|| denotes the Euclidean norm.
C 
C     W       is a DOUBLE PRECISION array of length N(2M+1)+2M used as
C             workspace for LBFGS. This array must not be altered by the
C             user.
C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      DOUBLE PRECISION :: GNORM, YY, SQ, YR, betalbfgs
      INTEGER CP, I, INMC, IYCN, ISCN
      INTEGER :: POINT, ISPT, IYPT, ITER, BOUND, NPT
      DOUBLE PRECISION :: STP, STP1, YS

      DOUBLE PRECISION, ALLOCATABLE :: XNEW(:),GNEW(:)
      DOUBLE PRECISION :: SLENGTH, OVERLAP, DOT1, DOT2, FNEW
      INTEGER :: NDECREASE

      character(len=512) :: tmpstring

      allocate(W(N*(2*M+1)+2*M), DIAG(N))
      allocate(XNEW(N),GNEW(N))

C
C
C     INITIALIZE
C     ----------
C

      CALL POT_AND_GRAD(X,F,G)

      ITER= 0
      POINT= 0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF(DIAGCO) THEN
C         DO I=1,N
C            IF (DIAG(I).LE.0.0D0) THEN
C              print *,'LBFGS: non-diagonal element is zero'
C              stop
C            END IF
C         END DO
C      ELSE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         DO I=1,N
            DIAG(I)= DGUESS
         END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
C     ---------------------------------------
C     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
C         OTHER TEMPORARY INFORMATION.
C     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
C     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
C         IN THE FORMULA THAT COMPUTES H*G.
C     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
C         STEPS.
C     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
C         GRADIENT DIFFERENCES.
C
C     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
C     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
C

      ISPT= N+2*M
      IYPT= ISPT+N*M     
      DO I=1,N
         W(ISPT+I)= -G(I)*DIAG(I)
         W(I)= -G(I)*DIAG(I)
      END DO
!      WRITE(*,*) 'DEBUG - Gradient vector:'
!      WRITE(*,*) G
      GNORM= DSQRT(DDOT(N,G,1,G,1))
      STP1= MIN(1.0D0/GNORM,GNORM)


      call opt_event_initialisation_done(f,gnorm,eps)

C
C    --------------------
C     MAIN ITERATION LOOP
C    --------------------
C

      main: DO
      ! Iteration is converged:
      IF (GNORM .LE. EPS) THEN
              call opt_event_opt_converged_successfully(iter,x,f,gnorm)
              exit
      END IF
      IF (ITER.EQ.MAXITER) THEN
              call opt_event_opt_max_iteration(iter,x,f,gnorm)
              exit 
      END IF

      IF (ITER .eq. 0) THEN
              call opt_event_first_iteration(x,f,gnorm,eps)
      end if

      ITER= ITER+1
      BOUND=ITER-1

      IF (ITER.EQ.1) THEN

         STP=STP1

      ELSE

         IF (ITER .GT. M) BOUND=M
         YS= DDOT(N,W(IYPT+NPT+1:),1,W(ISPT+NPT+1:),1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF(.NOT.DIAGCO) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         YY= DDOT(N,W(IYPT+NPT+1:),1,W(IYPT+NPT+1:),1)
         IF (YY.EQ.0.0D0) YY=1.0D0
         DO I=1,N
            DIAG(I)= YS/YY
         END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ELSE
C         input diagonal hessian elements here
C      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

*** set the diagonal elements here if diagco is true
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF(DIAGCO) THEN
C        DO I=1,N
C           IF (DIAG(I).LE.0.0D0) THEN
C              print *,'LBFGS: non-diagonal element is zero'
C              stop
C           END IF
C        END DO
C      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
C     "Updating quasi-Newton matrices with limited storage",
C     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
C     ---------------------------------------------------------
C
      CP= POINT
      IF (POINT.EQ.0) CP=M

      W(N+CP)= 1.0D0/YS

      DO I=1,N
         W(I)= -G(I)
      END DO

      CP= POINT

      DO I= 1,BOUND
         CP=CP-1
         IF (CP.EQ. -1) CP=M-1
         SQ= DDOT(N,W(ISPT+CP*N+1:),1,W,1)
         INMC=N+M+CP+1
         IYCN=IYPT+CP*N
         W(INMC)= W(N+CP+1)*SQ
         CALL DAXPY(N,-W(INMC),W(IYCN+1:),1,W,1)
      END DO

      DO I=1,N
         W(I)=DIAG(I)*W(I)
      END DO

      DO I=1,BOUND
         YR= DDOT(N,W(IYPT+CP*N+1:),1,W,1)
         betalbfgs= W(N+CP+1)*YR
         INMC=N+M+CP+1
         betalbfgs= W(INMC)-betalbfgs
         ISCN=ISPT+CP*N
         CALL DAXPY(N,betalbfgs,W(ISCN+1:),1,W,1)
         CP=CP+1
         IF (CP.EQ.M)CP=0
      END DO

      STP=1.0D0

      END IF



C     STORE THE NEW SEARCH DIRECTION
C     ------------------------------
       DO I=1,N
          W(ISPT+POINT*N+I)= W(I)
       END DO

****************************************************
***   INSTEAD OF LINE SEARCH
***   SET NEW STP, X, F, G
****************************************************
      DOT1=SQRT(DDOT(N,G,1,G,1))
      DOT2=SQRT(DDOT(N,W,1,W,1))
      OVERLAP=0.0D0
      IF (DOT1*DOT2 .NE. 0.0D0) THEN
         OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)
      END IF

      IF (OVERLAP .GT. 0.0D0) THEN
         DO I=1,N
            W(ISPT+POINT*N+I)=-W(I)
         END DO
      END IF

      DO I=1,N
         W(I)=G(I)
      END DO

      SLENGTH=0.0D0
      DO I=1,N
         SLENGTH=SLENGTH+W(ISPT+POINT*N+I)**2
      ENDDO
      SLENGTH=SQRT(SLENGTH)


      IF (STP*SLENGTH.GT.MAXMOVE) STP=MAXMOVE/SLENGTH


*** cycle to determine the step
      NDECREASE=0
      stepcycle: DO

      DO I=1,N
         XNEW(I)=X(I)+STP*W(ISPT+POINT*N+I)
      ENDDO

      CALL POT_AND_GRAD(XNEW,FNEW,GNEW)


      IF (FNEW.EQ.0.0D0) FNEW=1.0D-100 ! to prevent divide by zero
      IF ((FNEW-F)/ABS(FNEW).LE.MAXERISE) THEN
         X=XNEW
         F=FNEW
         G=GNEW
         EXIT
      ELSE
         IF (NDECREASE.GT.5) THEN
                 tmpstring = 'LBFGS: failed after 5 decreases'
                 call opt_event_iteration_error(ITER, tmpstring)
                 IFAIL=IFAIL+1
                 RETURN
         ELSE
            NDECREASE=NDECREASE+1
            STP=STP/10.0D0
         END IF
      END IF

      END DO stepcycle
****************************************************

C
C     COMPUTE THE NEW STEP AND GRADIENT CHANGE 
C     -----------------------------------------
C
      NPT=POINT*N

      DO I=1,N
         W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
         W(IYPT+NPT+I)= G(I)-W(I)
      END DO

      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0

C     TERMINATION TEST
C     ----------------

!      WRITE(*,*) 'DEBUG - Gradient vector:'
!      WRITE(*,*) G
      GNORM= DSQRT(DDOT(N,G,1,G,1))

      call opt_event_iteration_finished(iter,x,f,gnorm)
      END DO main

      END SUBROUTINE LBFGS
C
C     LAST LINE OF SUBROUTINE LBFGS
C
C

C   ----------------------------------------------------------
C
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision, dimension(:) :: dx, dy
      double precision :: da
      integer :: i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue

      end subroutine daxpy
C
C
C   ----------------------------------------------------------
C
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision, dimension(:) :: dx, dy
      double precision :: dtemp
      integer :: i, incx, incy, ix, iy, m, mp1, n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      END FUNCTION DDOT

****************************************************************************
      end module lbfgs_module

