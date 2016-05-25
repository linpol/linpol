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
Module opt_interface
use lbfgs_module, only: lbfgs_init => optim_init, lbfgs_opt => optimize
!use opt_events, only: set_linpol
implicit none
private
!
! M. Herbst:
! common interface to optimisation routines. Also provides events to be called
! by the optimiser and an enum to simplify selection.
! 
! TODO: Use a generic linpol data type

type, public :: opt
        integer,private :: opt_type
        contains
                procedure :: init => opt_init
                procedure :: opt_name => opt_opt_name
                procedure :: optimise => opt_optimise
end type

! define opt parameters:
integer, public, parameter :: OPT_NONE = 0
integer, public, parameter :: OPT_LBFGS = 1

public :: opt_type_from_string, opt_string_from_type

! NEW OPT METHODS:
! If a new optimisation method should be added you have to do the following:
!               - amend all routines in this module, especially the parameter
!               above and the pes_pot_name, pes_type_from_string methods (ADD_HERE)
!
! ===================================================
! #- Current optimisations -#
! ===========================
!
! LBFGS:
!       some properies
!
! ==================================================

contains
subroutine opt_init(this, opt_type)
use phys_const
        class(opt),intent(out) :: this
        integer, intent(in) :: opt_type

        this%opt_type = opt_type

        SELECT CASE(opt_type)
                case(OPT_LBFGS)
                        return !nothing to do
                !ADD HERE
                case default
                        stop "Cannot initialise an empty pes type"
        end select
end subroutine

subroutine opt_optimise(this,x,f,g,msave,eps,maxiterations)
        class(opt),intent(in) :: this
        double precision, dimension(:),intent(in) :: x !inital coords -- TODO: replace by linpol later

        double precision, dimension(:),intent(out) :: g !final gradient
        double precision, intent(out) :: f !final function value


        integer,intent(in) :: msave ! Number of saved steps, between 3 and 7
        double precision, intent(in) :: eps ! Epsilon for stopping criterion (if grad is under it is stopped)
        integer, intent(in) :: maxiterations ! Max number of iterations until minimisation aborts

        
        ! set_linpol(linpol)


        SELECT CASE(this%opt_type)
                case(OPT_LBFGS)
                        call lbfgs_init( 4,msave,eps,maxiterations ) !reset lbfgs module for next iteration
                        !choose to call these functions
                        CALL lbfgs_opt ( x, f, g, msave,eps,maxiterations)
                !ADD HERE
        end select
end subroutine

subroutine opt_opt_name(this,opt_string)       ! get the name of pot object
        class(opt), intent(in) :: this                   !selected pes
        character(len=8), intent(out) :: opt_string
        call opt_string_from_type(this%opt_type,opt_string)
end subroutine

! ==================================================

subroutine opt_type_from_string(opt_string,opt_type)    ! get the pot object integer from string
        character(len=8), intent(in) :: opt_string
        integer,intent(out) :: opt_type

        select case (opt_string)
                case("LBFGS")
                        opt_type = OPT_LBFGS
                !ADD HERE
                case default
                        opt_type = OPT_NONE
        end select
end subroutine

subroutine opt_string_from_type(opt_type,opt_string)    ! get the pot name from pot object integer
        character(len=8), intent(out) :: opt_string
        integer,intent(in) :: opt_type

        SELECT CASE(opt_type)
                case(OPT_LBFGS)
                        opt_string = "LBFGS"
                !ADD HERE
                case default
                        opt_string = "NONE"
        end select
end subroutine

end module
