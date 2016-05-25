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
module opt_events
use pots
use grads
use parameters, only: root, pes_obj
implicit none
private
!
! M. Herbst:
! events for minimisation procedures.

!public :: set_linpol
public :: opt_event_next_pot_grad, opt_event_opt_restarted, opt_event_opt_failed, opt_event_initialisation_done
public :: opt_event_first_iteration, opt_event_iteration_error, opt_event_iteration_finished
public :: opt_event_opt_converged_successfully, opt_event_opt_max_iteration

double precision, save :: lowestF, lowestGNORM

!type(linpolobj) :: tmplinpol !TODO: save the linpol object here

contains
!subroutine set_linpol(linpol)
!       tmplinpol = linpol
!end subroutine

subroutine opt_event_next_pot_grad(x,f,g)
        double precision, dimension(:), intent(in) :: x
        double precision, intent(out) :: f
        double precision, dimension(:), intent(out) :: g

        if (singleBeadMode) then
                !single bead:
                call pes_obj%pot_grad(x,gradeps,f,g)
        else
                !multiple bead:
                CALL totpotx ( x, f )
                CALL totgradx ( x, g )
        end if
end subroutine

!=========================================================

subroutine opt_event_opt_restarted(failcount)
        use mpi
        integer,intent(in) :: failcount !how many iterations failed
        integer :: myrank,error

        CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)
        if (myrank .ne. root) return

        print *,'LBFGS: lbfgs failed once and is restarted'
end subroutine

subroutine opt_event_opt_failed(failcount)
        use mpi
        integer,intent(in) :: failcount !how many iterations failed
        integer :: myrank,error

        CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)
        if (myrank .ne. root) return

        print *,'LBFGS: lbfgs failed again and is terminated'
end subroutine

!=========================================================

subroutine opt_event_initialisation_done(f,gnorm,eps)
        use mpi
        double precision, intent(in) :: f
        double precision, intent(in) :: gnorm 
        double precision, intent(in) :: eps
        integer :: myrank,error

        CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)
        if (myrank .ne. root) return

        ! To report on initial geometry fed into L-BFGS
        WRITE(99,*)
        WRITE(*,*)
        WRITE(*,*) 'Initial GNORM: ', GNORM, ' hartree/bohr'
        WRITE(*,*) 'Initial energy: ', F, ' hartree'
        WRITE(99,*) 'Initial GNORM: ', GNORM, ' hartree/bohr'
        WRITE(99,*) 'Initial energy: ', F, ' hartree'
end subroutine

subroutine opt_event_first_iteration(x,f,gnorm,eps)
        use mpi
        use parameters,only: latest_gnorm
        use parameters,only: savefile, lowestgradfile, lowestfuncfile
        double precision, intent(in) :: f
        double precision, dimension(:), intent(in) :: x
        double precision, intent(in) :: gnorm 
        double precision, intent(in) :: eps
        integer :: myrank,error
        
        latest_gnorm = GNORM         !TODO MFH: Don't really like this:
        CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)
        if (myrank .ne. root) return

        ! Better have data structure containing supporting iteration
        ! information, which should also save this 


        WRITE(99,*)
        WRITE(99,'(A)') 'Guide to letters at end of some of the iteration lines:'
        WRITE(99,'(A,A)') ' s: geometry saved to ', savefile
        WRITE(99,'(A,A)') ' g: geometry with lowest gradient saved to ', lowestgradfile
        WRITE(99,'(A,A)') ' f: geometry with lowest function value saved to ', lowestfuncfile
        WRITE(99,*)
        WRITE(*,*)
end subroutine

subroutine opt_event_iteration_error(iter,message)
        use mpi
        integer,intent(in) :: iter
        character(len=512),intent(in) :: message
        integer :: myrank,error

        CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)
        if (myrank .ne. root) return

        write(99,*) trim(message)
        write(99,*) 'NUMBER OF ITERATIONS IS',ITER
end subroutine

subroutine opt_event_iteration_finished(iter,x,f,gnorm)
        use mpi
        use parameters, only: writeline, savestate

        !TODO: temporary:
        USE funcs, only: onetothree, xyzout
        USE prep, only: d, p, beta, atomtags
        !//// temporary

        integer,intent(in) :: iter
        double precision, dimension(:), intent(in) :: x
        double precision, intent(in) :: f
        double precision, intent(in) :: gnorm 
        integer :: myrank,error

        CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)
        if (myrank .ne. root) return

        if ( MOD(iter,writeline) == 0 ) then
                WRITE(99,2222) iter, GNORM, F
2222            FORMAT ('Iter:',I5,'  GNORM:',ES12.5,'  ENGY:',F0.12,' ',$)
                WRITE(*,2221) iter, GNORM, F
2221            FORMAT ('\rIter:', I5, '  GNORM:', ES12.5, '  ENGY:', F0.8, $)
                ! Dollar sign at end of previous command is included to suppress the automatic new line
                ! \r is carriadge return, ie goes back to beginning of line
        end if

        IF (iter == 1) THEN
                lowestF = F
                lowestGNORM = GNORM
        END IF
        
        if ( mod(iter,savestate) == 0 ) then
                CALL xyzout(x, size(x)/(d*p) ,d,p,beta,atomtags,savefile)
                WRITE(99,2223)
2223            FORMAT ('s', $ )
        end if

        if ( GNORM .LT. lowestGNORM ) THEN
                CALL xyzout(x, size(x)/(d*p) ,d,p,beta,atomtags,lowestgradfile)
                WRITE(99,2224)
2224            FORMAT ('g', $ )
                lowestGNORM = GNORM
        END IF

        IF ( F .LT. lowestF ) THEN
                CALL xyzout(x, size(x)/(d*p) ,d,p,beta,atomtags,lowestfuncfile)
                WRITE(99,2225)
2225            FORMAT ('f', $ )
                lowestF = F
        END IF
        WRITE(99,*)
end subroutine

subroutine opt_event_opt_converged_successfully(iter,x,f,gnorm)
        use mpi
        use parameters,only: latest_gnorm
        integer,intent(in) :: iter
        double precision, dimension(:), intent(in) :: x
        double precision, intent(in) :: f
        double precision, intent(in) :: gnorm 
        integer :: myrank,error

        latest_gnorm = gnorm
        CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)
        if (myrank .ne. root) return

        WRITE(*,*)
        WRITE(*,*)
        WRITE(*,'(A)') '*************************************'
        WRITE(*,'(A)') '*** L-BFGS CONVERGED SUCCESSFULLY ***'
        WRITE(*,'(A)') '*************************************'
        WRITE(*,*)
        WRITE(99,*)
        WRITE(99,'(A)') '*************************************'
        WRITE(99,'(A)') '*** L-BFGS CONVERGED SUCCESSFULLY ***'
        WRITE(99,'(A)') '*************************************'
        WRITE(99,*)
        WRITE(*,'(A,I5)') 'Number of iterations was: ', ITER
        WRITE(99,'(A,I5)') 'Number of iterations was: ', ITER
        WRITE(*,*) 'Closing GNORM was: ', GNORM, ' hartree/bohr'
        WRITE(99,*) 'Closing GNORM was: ', GNORM, ' hartree/bohr'
        WRITE(*,*) 'Closing energy was: ', F, ' hartree'
        WRITE(99,*) 'Closing energy was: ', F, ' hartree'
        WRITE(*,*)
        WRITE(99,*)
end subroutine

subroutine opt_event_opt_max_iteration(iter,x,f,gnorm)
        use mpi
        use parameters,only: latest_gnorm
        integer,intent(in) :: iter
        double precision, dimension(:), intent(in) :: x
        double precision, intent(in) :: f
        double precision, intent(in) :: gnorm 
        integer :: myrank,error

        latest_gnorm = gnorm
        CALL MPI_Comm_rank ( MPI_COMM_WORLD, myrank, error)
        if (myrank .ne. root) return

        WRITE(*,*)
        WRITE(*,*)
        WRITE(*,'(A)') '*********************************************'
        WRITE(*,'(A)') '*** MAX ITERATIONS REACHED: NOT CONVERGED ***'
        WRITE(*,'(A)') '*********************************************'
        WRITE(*,*)
        WRITE(99,*)
        WRITE(99,'(A)') '*********************************************'
        WRITE(99,'(A)') '*** MAX ITERATIONS REACHED: NOT CONVERGED ***'
        WRITE(99,'(A)') '*********************************************'
        WRITE(99,*)
        WRITE(*,'(A,I5)') 'Number of iterations was: ', ITER
        WRITE(99,'(A,I5)') 'Number of iterations was: ', ITER
        WRITE(*,*) 'Closing GNORM was: ', GNORM, ' hartree/bohr'
        WRITE(99,*) 'Closing GNORM was: ', GNORM, ' hartree/bohr'
        WRITE(*,*) 'Closing energy was: ', F, ' hartree'
        WRITE(99,*) 'Closing energy was: ', F, ' hartree'
        WRITE(*,*)
        WRITE(99,*)
end subroutine

end module
