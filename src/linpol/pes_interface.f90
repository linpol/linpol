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
Module pes_interface
!use pes_shell, only: whbb_init => pes_init, whbb_pot => f, whbb_grad => grad, whbb_hess => fhess2
use mod_ttm3f
implicit none
private
!
! M. Herbst:
! common interface to the various potentials. Provides initialiser and
! enum to simplify selection.
! NOTE: Assumes that instanton is saved in a 1-D fortran array of size 9*Nw
! The actual potential, gradients and hessian can only be calculated properly,
! if the coordinate input x uses the correct convention of the potential itself.
! Hence before getting the potential for the first time the coordinates should
! be transformed into the convention used by the pes potential using the vec_to_pes,
! mat_to_pes subroutines. Before doing any output and also after the
! minimisation / other routines using the potential the transformation should be
! reversed using pes_mat_from_pes, pes_vec_from_pes.

type, public :: pes
        integer,private :: pes_type
        integer,private :: Nw
        double precision,private :: EnergyUnit = 1.d0 !1Hartree in units of energy used by the pes
        double precision,private :: LengthUnit = 1.d0 !1bohr in units of length used by the pes
        contains
                procedure :: init => pes_init
                procedure :: pot => pes_pot
                procedure :: grad => pes_grad
                procedure :: pot_grad => pes_pot_grad
                procedure :: hess => pes_hess
                procedure :: pes_name => pes_pes_name
                procedure :: vec_from_pes => pes_vec_from_pes
                procedure :: vec_to_pes => pes_vec_to_pes
                procedure :: atomvec_from_pes => pes_atomvec_from_pes
                procedure :: atomvec_to_pes => pes_atomvec_to_pes
                procedure :: mat_from_pes => pes_mat_from_pes
                procedure :: mat_to_pes => pes_mat_to_pes
end type

! define pes parameters:
integer, public, parameter :: PES_NONE = 0
integer, public, parameter :: PES_TTM3_F = 1
integer, public, parameter :: PES_MBPOL = 2

public :: pes_type_from_string, pes_string_from_type

! NEW POTENTIALS:
! If a new potential should be added you have to do the following:
!               - amend all routines in this module, especially the parameter
!               above and the pes_pot_name, pes_type_from_string methods (ADD_HERE)
!               - unless the new potential you are adding uses the same convention
!               for ordering the atoms and has the same units as HBB2 and WHBB, you
!               will need to add cases carefully to all of the subroutines below
!               - create a folder of the chosen potential name to data

! NOTE: this program works in atomic units => all return values need to be
! converged using the physical constants!

! ===================================================
! #- Current potentials -#
! ========================

! MBPOL:
! energy is in kcal/mol
! length is in A
! gradient is in kcal/mol/A
! Coordinates stored as:
!       OA
!       H1              H1 - OA - H2
!       H2
!       OB
!       H3              H3 - OB - H4
!       H4

! TTM3_F:
! energy is in kcal/mol
! length is in A
! gradient is in kcal/mol/A
! Coordinates stored as
!       OA
!       OB              H1 - OA - H2
!       H1
!       H2
!       H3              H3 - OB - H4
!       H4

! ==================================================
contains
subroutine pes_init(this, pes_type, Nw)
use phys_const
        class(pes),intent(out) :: this
        integer, intent(in) :: pes_type
        integer, intent(in) :: Nw       !number of water molecules

        this%Nw = Nw
        this%pes_type = pes_type

        SELECT CASE(pes_type)
                case(PES_TTM3_F)
                        call init_ttm3f(Nw)
                        this%EnergyUnit = aukcal
                        this%LengthUnit = auang
                case(PES_MBPOL)
                        this%EnergyUnit = aukcal
                        this%LengthUnit = auang
                !ADD HERE
                case default
                        stop "Cannot initialise an empty pes type"
        end select
end subroutine

subroutine pes_vec_from_pes(this, v)
        ! transform a vector (like Position) from coordinate convention used by
        ! the pes potential to coordinate convention used by the program
        class(pes), intent(in) :: this                          ! selected pes
        double precision, dimension(9*this%Nw), intent(inout) :: v
        double precision, allocatable,dimension(:) :: tmp
        integer water
        integer Hcoords
        integer Ocoords

        SELECT CASE(this%pes_type)
                case(PES_TTM3_F)
                ! INPUT ORDER (TTM3-F CONVENTION):
                !       OA
                !       OB              H1 - OA - H2
                !       H1
                !       H2
                !       H3              H3 - OB - H4
                !       H4
                        allocate( tmp(3*this%Nw) )   !temporary storage for oxygens
                        tmp = v(1:3*this%Nw)
                        v(1:6*this%Nw) = v(3*this%Nw+1:9*this%Nw)       !pop Hs to front
                        v(6*this%Nw+1:9*this%Nw) = tmp
                        deallocate( tmp )
                ! OUTPUT ORDER (PROGRAM CONVENTION):
                !       H1
                !       H2              H1 - OA - H2
                !       H3
                !       H4
                !       OA              H3 - OB - H4
                !       OB
                case(PES_MBPOL)
                ! INPUT ORDER (MBPOL CONVENTION):
                !       OA
                !       H1              H1 - OA - H2
                !       H2
                !       OB
                !       H3              H3 - OB - H4
                !       H4
                    allocate( tmp(9*this%Nw) )   !temporary storage for coordinates
                    do water=1,this%Nw	
                      do Ocoords=1,3
                        tmp( (this%Nw*6)+(3*(water-1))+Ocoords ) = v( (9*(water-1))+Ocoords ) ! move the oxygens to the end
                      end do
                      do Hcoords=1,6
                        tmp( (6*(water-1))+Hcoords ) = v( 3+(9*(water-1))+Hcoords ) ! move the hydrogens in front of the oxygens
                      end do
                    end do
                    v = tmp
                    deallocate( tmp )
                ! OUTPUT ORDER (PROGRAM CONVENTION):
                !       H1
                !       H2              H1 - OA - H2
                !       H3
                !       H4
                !       OA              H3 - OB - H4
                !       OB

                !ADD HERE
        end select
end subroutine

subroutine pes_vec_to_pes(this, v)
        ! transform a 9*Nw vector (like Position) from coordinate convention used by
        ! the program to coordinate convention used by the pes potential
        class(pes), intent(in) :: this                          ! selected pes
        double precision, dimension(9*this%Nw), intent(inout) :: v
        double precision, allocatable,dimension(:) :: tmp
        integer water
        integer Hcoords
        integer Ocoords

        SELECT CASE(this%pes_type)
                case(PES_TTM3_F)
                        allocate( tmp(3*this%Nw) )   !temporary storage for oxygens
                        tmp = v(6*this%Nw+1:9*this%Nw)
                        v(3*this%Nw+1:9*this%Nw) = v(1:6*this%Nw)       !push back hydrogens
                        v(1:3*this%Nw) = tmp                            !oxygens at front
                        deallocate( tmp )
                case(PES_MBPOL)
                ! MBPOL: INPUT ORDER (PROGRAM CONVENTION):
                !       H1
                !       H2              H1 - OA - H2
                !       H3
                !       H4
                !       OA              H3 - OB - H4
                !       OB
                    allocate( tmp(9*this%Nw) )   !temporary storage for coordinates
                    do water=1,this%Nw	
                      do Ocoords=1,3
                        tmp( (9*(water-1))+Ocoords ) = v( (this%Nw*6)+(3*(water-1))+Ocoords ) ! move the oxygens from the end
                      end do
                      do Hcoords=1,6
                        tmp( 3+(9*(water-1))+Hcoords ) = v( (6*(water-1))+Hcoords ) ! move the hydrogens from in front of the oxygens
                      end do
                    end do
                    v = tmp
                    deallocate( tmp )
                ! MBPOL: OUTPUT ORDER (PES CONVENTION):
                !       OA
                !       H1              H1 - OA - H2
                !       H2
                !       OB
                !       H3              H3 - OB - H4
                !       H4

                !ADD HERE
        end select
end subroutine

subroutine pes_atomvec_from_pes(this, v)
        ! transform a 3*Nw vector (like atom masses) from coordinate convention used by
        ! the pes potential to coordinate convention used by the program
        class(pes), intent(in) :: this                          ! selected pes
        double precision, dimension(3*this%Nw), intent(inout) :: v
        double precision, allocatable,dimension(:) :: tmp
        integer water
        integer Hatoms

        SELECT CASE(this%pes_type)
                case(PES_TTM3_F)
                        allocate( tmp(this%Nw) )   !temporary storage for oxygens
                        tmp = v(1:this%Nw)
                        v(1:2*this%Nw) = v(this%Nw+1:3*this%Nw)       !pop Hs to front
                        v(2*this%Nw+1:3*this%Nw) = tmp
                        deallocate( tmp )
                case(PES_MBPOL)
                ! MBPOL: INPUT ORDER (PES CONVENTION):
                !       OA
                !       H1              H1 - OA - H2
                !       H2
                !       OB
                !       H3              H3 - OB - H4
                !       H4
                        allocate( tmp(3*this%Nw) )   !temporary storage for coordinates
                        do water=1,this%Nw
                          tmp( (this%Nw*2)+water ) = v( 3*(water-1)+1 ) ! move the oxygens to the end
                          do Hatoms=1,2
                            tmp( (2*(water-1))+Hatoms ) = v( 1+(3*(water-1))+Hatoms ) ! move the hydrogens in front of the oxygens
                          end do
                        end do
                        v = tmp
                        deallocate( tmp )
                ! MBPOL: OUTPUT ORDER (PROGRAM CONVENTION):
                !       H1
                !       H2              H1 - OA - H2
                !       H3
                !       H4
                !       OA              H3 - OB - H4
                !       OB

                !ADD HERE
        end select
end subroutine

subroutine pes_atomvec_to_pes(this, v)
        ! transform a 3*Nw vector (like atom masses) from coordinate convention used by
        ! the program to coordinate convention used by the pes potential
        class(pes), intent(in) :: this                          ! selected pes
        double precision, dimension(3*this%Nw), intent(inout) :: v
        double precision, allocatable,dimension(:) :: tmp
        integer water
        integer Hatoms

        SELECT CASE(this%pes_type)
                case(PES_TTM3_F)
                        allocate( tmp(this%Nw) )   !temporary storage for oxygens
                        tmp = v(2*this%Nw+1:3*this%Nw)
                        v(this%Nw+1:3*this%Nw) = v(1:2*this%Nw)       !push back hydrogens
                        v(1:this%Nw) = tmp                            !oxygens at front
                        deallocate( tmp )
                case(PES_MBPOL)
                ! MBPOL: INPUT ORDER (PROGRAM CONVENTION):
                !       H1
                !       H2              H1 - OA - H2
                !       H3
                !       H4
                !       OA              H3 - OB - H4
                !       OB
                        allocate( tmp(3*this%Nw) )   !temporary storage for coordinates
                        do water=1,this%Nw
                          tmp( 3*(water-1)+1 ) = v( (this%Nw*2)+water ) ! move the oxygens from the end
                          do Hatoms=1,2
                            tmp( 1+(3*(water-1))+Hatoms ) = v( (2*(water-1))+Hatoms ) ! move the hydrogens from in front of the oxygens
                          end do
                        end do
                        v = tmp
                        deallocate( tmp )
                ! MBPOL: OUTPUT ORDER (PES CONVENTION):
                !       OA
                !       H1              H1 - OA - H2
                !       H2
                !       OB
                !       H3              H3 - OB - H4
                !       H4
                
                !ADD HERE
        end select
end subroutine

subroutine pes_mat_from_pes(this, m)
        ! transform a matrix (like Hessian) from coordinate convention used by
        ! the pes potential to coordinate convention used by the program
        class(pes), intent(in) :: this                          ! selected pes
        double precision, dimension(9*this%Nw,9*this%Nw), intent(inout) :: m
        double precision, allocatable,dimension(:,:) :: tmp, tmp2
        integer :: water, Ocoords, Hcoords, i, j, dim
        integer, dimension(9*this%Nw) :: postranslate

        SELECT CASE(this%pes_type)
                case(PES_TTM3_F)
                        ! We have blocks:       We want blocks like
                        !     ( O:O, O:H )          ( H:H, H:O )
                        !     ( H:O, H:H )          ( O:H, O:O )
                        !
                        ! hence we temporarily save the row block ( O:O, O:H )
                        allocate( tmp(3*this%Nw,9*this%Nw) )
                        tmp = m(1:3*this%Nw,1:9*this%Nw)
                        !
                        ! and the block H:O
                        allocate( tmp2(6*this%Nw,3*this%Nw) )
                        tmp2 = m(3*this%Nw+1:9*this%Nw,1:3*this%Nw)
                        !
                        ! and permute H:H to the upper left
                        m(          1:6*this%Nw,          1:6*this%Nw) = m(3*this%Nw+1:9*this%Nw,3*this%Nw+1:9*this%Nw)
                        !
                        ! next we place H:O right to it
                        m(          1:6*this%Nw,6*this%Nw+1:9*this%Nw) = tmp2(       1:6*this%Nw,         1:3*this%Nw)
                        !
                        ! finally we put the (O:O) and (O:H) blocks down where they belong:
                        m(6*this%Nw+1:9*this%Nw,6*this%Nw+1:9*this%Nw) = tmp(        1:3*this%Nw,          1:3*this%Nw)
                        m(6*this%Nw+1:9*this%Nw,          1:6*this%Nw) = tmp(        1:3*this%Nw,3*this%Nw+1:9*this%Nw)
                        !
                        ! and free the memory:
                        deallocate( tmp,tmp2)
                case(PES_MBPOL)
                        ! Create array that translates atomic array positions from pes convention to program convention
                        do water=1,this%Nw	
                          do Ocoords=1,3
                            postranslate((9*(water-1))+Ocoords) = (this%Nw*6)+(3*(water-1))+Ocoords
                          end do
                          do Hcoords=1,6
                            postranslate(3+(9*(water-1))+Hcoords) = (6*(water-1))+Hcoords
                          end do
                        end do
                        ! Allocate array to hold the new second derivatives
                        allocate( tmp(9*this%Nw,9*this%Nw) )
                        dim=9*this%Nw
                          do i=1,dim
                            do j=1,dim
                              tmp( postranslate(i),postranslate(j) ) = m(i,j)
                            end do
                          end do
                        m = tmp
                        deallocate( tmp )
                !ADD HERE
        end select
end subroutine

subroutine pes_mat_to_pes(this, m)
        ! transform a matrix (like Hessian) from coordinate convention used by
        ! the program to coordinate convention used by the pes potential
        class(pes), intent(in) :: this                          ! selected pes
        double precision, dimension(9*this%Nw,9*this%Nw), intent(inout) :: m
        double precision, allocatable,dimension(:,:) :: tmp,tmp2

        SELECT CASE(this%pes_type)
                case(PES_TTM3_F)
                        ! We have blocks:       We want blocks like
                        !     ( H:H, H:O )          ( O:O, O:H )
                        !     ( O:H, O:O )          ( H:O, H:H )
                        !
                        ! hence we temporarily save the row block ( O:H, O:O )
                        allocate( tmp(3*this%Nw,9*this%Nw) )
                        tmp = m(6*this%Nw+1:9*this%Nw,1:9*this%Nw)
                        !
                        ! and the block H:O
                        allocate( tmp2(6*this%Nw,3*this%Nw) )
                        tmp2 = m(1:6*this%Nw,6*this%Nw+1:9*this%Nw)
                        !
                        ! and permute H:H to the lower right
                        m(3*this%Nw+1:9*this%Nw,3*this%Nw+1:9*this%Nw) = m(1:6*this%Nw, 1:6*this%Nw)
                        !
                        ! next we place H:O left to it
                        m(3*this%Nw+1:9*this%Nw,1:3*this%Nw) = tmp2(1:6*this%Nw,1:3*this%Nw)
                        !
                        ! finally we put the rest down where it belongs:
                        m(1:3*this%Nw,1:3*this%Nw) = tmp(1:3*this%Nw,6*this%Nw+1:9*this%Nw)
                        m(1:3*this%Nw,3*this%Nw+1:9*this%Nw) = tmp(1:3*this%Nw,1:6*this%Nw)
                        !
                        ! and free the memory:
                        deallocate( tmp,tmp2)
                case(PES_MBPOL)
                        write(*,*) ''
                        write(*,*) '***********************************************'
                        write(*,*) 'pes_mat_to_pes not yet implemented for MB-pol!'
                        write(*,*) '***********************************************'
                        write(*,*) ''
                        stop 1
                !ADD HERE
        end select
end subroutine

subroutine pes_pot(this, x, pot)    !calculate potential
        !NOTE: We assume the correct coordinate convention is used and input is
        !in atomic units!
        class(pes), intent(in) :: this                          ! selected pes
        double precision, dimension(9*this%Nw), intent(in) :: x	! In units of Bohr
        double precision, intent(out) :: pot                    ! potential
        double precision, dimension(:),allocatable :: grad ! to be discarded, neccessary for ttm3f
        double precision, dimension(9*this%Nw) :: xunit

        !Forward transform units (from Bohr to potential's native units):
        xunit = x*this%LengthUnit

        SELECT CASE(this%pes_type)
                case(PES_TTM3_F)
                        allocate(grad(9*this%Nw)) !TODO: Get rid of this somehow!!
                        call ttm3f(this%Nw,xunit,grad,pot, .false.)    !the false disables gradient calculation
                        deallocate(grad)
                case(PES_MBPOL)
                        call calcpot(this%Nw,pot,xunit)
                        !Potential is returned in units of kcal/mol; these are converted below, to Ha
                !ADD HERE
        end select

        !Backward transform units to Ha:
        pot=pot/this%EnergyUnit

end subroutine

subroutine pes_grad(this, x, eps, grad)    !calculate gradient
        !NOTE: We assume the correct coordinate convention is used and input is
        !in atomic units!
        ! NOTE ON UNITS:
        ! - This subroutine takes x in Bohr (the standard for the Linpol program) and then
        !   scales the geometry to the correct units for the potential call (e.g. A, in the
        !   case of HBB2pol).  There is thus no need for the user to convert the geometry
        !   manually here when adding a new potential.
        ! - This subroutine takes eps in Bohr.  If the potential being used required eps to
        !   be in different units, the user will need to make this change when adding the new
        !   potential.
        class(pes), intent(in) :: this                            ! selected pes
        double precision,dimension(9*this%Nw),intent(in) :: x     ! coordinates, in Bohr
        double precision, intent(in) :: eps                       ! epsilon from user (au)
        double precision,dimension(9*this%Nw),intent(out) :: grad ! gradient
        double precision :: pot                      		  ! to be discarded, necessary for ttm3f, HBB2pol and MBpol
        double precision, dimension(9*this%Nw) :: xunit		  ! in potential's native units
        double precision :: GNORM

        !Forward transform units (from Bohr to potential's native units):
        xunit = x*this%LengthUnit

        SELECT CASE(this%pes_type)
                case(PES_TTM3_F)
                        call ttm3f(this%Nw,xunit,grad,pot)
                        !Note: pot discarded
                case(PES_MBPOL)
                        call calcpotg(this%Nw,pot,xunit,grad)
                        !Note: pot discarded
                        !Note: no eps option, as gradient is exact
                        !Gradient is returned in units of kcal/mol/A; these are converted below
                !ADD HERE
        end select

        !Backward transform units to Ha/bohr:
        grad = (grad/this%EnergyUnit)*this%LengthUnit

end subroutine

subroutine pes_pot_grad(this, x, eps, pot, grad)    !calculate gradient
        !NOTE: We assume the correct coordinate convention is used and input is
        !in atomic units!
        ! use this in preference over separately pot and grad, since quicker!
        class(pes), intent(in) :: this                            ! selected pes
        double precision,dimension(9*this%Nw),intent(in) :: x     ! coordinates, in Bohr
        double precision, intent(in) :: eps                       ! epsilon from user (au)
        double precision,dimension(9*this%Nw),intent(out) :: grad ! gradient
        double precision, intent(out) :: pot                      ! potential
        double precision, dimension(9*this%Nw) :: xunit		  ! in potential's native units
        double precision :: GNORM
        
        !Forward transform units (from Bohr to potential's native units):
        xunit = x*this%LengthUnit

        SELECT CASE(this%pes_type)
                case(PES_TTM3_F)
                        call ttm3f(this%Nw,xunit,grad,pot)
                case(PES_MBPOL)
                        call calcpotg(this%Nw,pot,xunit,grad)
                        !Note: no eps option, as gradient is exact
                        !Potential is returned in units of kcal/mol; these are converted below, to Ha
                        !Gradient is returned in units of kcal/mol/A; these are converted below, to Ha/bohr
                !ADD HERE
        end select

        !Backward transform units to Ha/bohr and Ha, respectively:
        grad = grad/this%EnergyUnit*this%LengthUnit
        pot=pot/this%EnergyUnit

end subroutine

subroutine pes_five_point_stencil_gradient_from_potential(this, x, eps, grad)
        ! M. Herbst:
        ! use the five-point stencil formula to calculate the gradient
        ! is correct up to fourth order
        ! NOTE ON UNITS:
        ! - This subroutine takes x in Bohr (the standard for the Linpol program) and feeds it
        !   into the pot function.  The pot function then deals with scaling the geometry
        !   to the correct units for the potential call (e.g. A, in the case of HBB2pol).  There
        !   is thus no need for the user to convert the geometry manually here.
        ! - This subroutine takes eps in Bohr and then uses it to generate slightly amended
        !   geometries, also in Bohr.  It then feeds these geometries into the pot function.
        !   The pot function then deals with scaling the geometries to the correct units for
        !   the potential call (e.g. A, in the case of HBB2pol).  There is thus no need for the
        !   user to convert the eps value manually.
        class(pes), intent(in) :: this                                      ! selected pes
        double precision,dimension(9*this%Nw),intent(in) :: x               ! coordinates
        double precision, intent(in) :: eps                                 ! epsilon for gradient calc
        double precision,dimension(9*this%Nw),intent(out) :: grad           ! gradient vector
        integer::dim,i
        double precision,dimension(9*this%Nw) :: xt
        double precision :: fb, fc, fd, fe
        !
        dim=size(x)
        do i=1,dim
          xt=x; xt(i)=xt(i)+2.0d0*eps
          call this%pot(xt,fb)
          xt=x; xt(i)=xt(i)+eps
          call this%pot(xt,fc)
          xt=x; xt(i)=xt(i)-eps
          call this%pot(xt,fd)
          xt=x; xt(i)=xt(i)-2.0d0*eps
          call this%pot(xt,fe)
          grad(i) = (-fb+8.0d0*fc-8.0d0*fd+fe)/(12.0d0*eps)
        end do
end subroutine

subroutine pes_standard_hessian_from_potential2(this, x, eps, hess)
        ! M. Herbst:
        ! use the standard formula to calculate a Hessian
        ! is correct up to second order only!
        ! Assumes all input is in atomic units
        class(pes), intent(in) :: this                                      ! selected pes
        double precision,dimension(9*this%Nw),intent(in) :: x               ! coordinates
        double precision, intent(in) :: eps                                 ! epsilon for hessian calc
        double precision,dimension(9*this%Nw,9*this%Nw),intent(out) :: hess ! hessian
        !
        integer i,j
        double precision,dimension(9*this%Nw) :: xt
        double precision :: fa, fb, fc, fd, fe, ff, fg
        !
        call this%pot(x,fa)
        do i=1,9*this%Nw
                xt=x; xt(i)=xt(i)+eps
                call this%pot(xt,fb)
                !
                xt=x; xt(i)=xt(i)-eps
                call this%pot(xt,fc)
                !
                hess(i,i) = (fb-2.0d0*fa+fc)/eps**2
                do j=1,i-1
                        xt=x; xt(i)=xt(i)+eps; xt(j)=xt(j)+eps
                        call this%pot(xt,fd)
                        !
                        xt=x; xt(i)=xt(i)-eps; xt(j)=xt(j)-eps
                        call this%pot(xt,fe)
                        !
                        xt=x; xt(j)=xt(j)+eps
                        call this%pot(xt,ff)
                        !
                        xt=x; xt(j)=xt(j)-eps
                        call this%pot(xt,fg)
                        !
                        ! Divide by epsilon in correct units!
                        hess(i,j)=0.5d0*(fd+fe-fb-fc-ff+2.0d0*fa-fg)/(eps)**2
                        hess(j,i)=hess(i,j)
                end do
        end do
end subroutine

subroutine pes_five_point_stencil_hessian_from_potential(this, x, eps, hess)
        ! M. Herbst:
        ! use the five-point stencil formula to calculate a Hessian
        ! based on fhess2, the 'alternative numerical hessian of the water potential'
        ! which was written by JOR in the original Linpol code
        ! is correct up to fourth order
        ! NOTE ON UNITS:
        ! - This subroutine takes x in Bohr (the standard for the Linpol program) and feeds it
        !   into the pot function.  The pot function then deals with scaling the geometry
        !   to the correct units for the potential call (e.g. A, in the case of HBB2pol).  There
        !   is thus no need for the user to convert the geometry manually here.
        ! - This subroutine takes eps in Bohr and then uses it to generate slightly amended
        !   geometries, also in Bohr.  It then feeds these geometries into the pot function.
        !   The pot function then deals with scaling the geometries to the correct units for
        !   the potential call (e.g. A, in the case of HBB2pol).  There is thus no need for the
        !   user to convert the eps value manually.
        class(pes), intent(in) :: this                                      ! selected pes
        double precision,dimension(9*this%Nw),intent(in) :: x               ! coordinates
        double precision, intent(in) :: eps                                 ! epsilon from user (au)
        double precision,dimension(9*this%Nw,9*this%Nw),intent(out) :: hess ! hessian
        integer::dim,i,j
        double precision,dimension(9*this%Nw) :: xt
        double precision :: fa, fb, fc, fd, fe
        !
        dim=size(x)
        call this%pot(x,fa)
        do i=1,dim
          xt=x; xt(i)=xt(i)+2.0d0*eps
          call this%pot(xt,fb)
          xt=x; xt(i)=xt(i)+eps
          call this%pot(xt,fc)
          xt=x; xt(i)=xt(i)-eps
          call this%pot(xt,fd)
          xt=x; xt(i)=xt(i)-2.0d0*eps
          call this%pot(xt,fe)
          hess(i,i)=(-fb+16.0d0*fc-30.0d0*fa+16.0d0*fd-fe)/(12.0d0*eps**2)
          do j=1,i-1
            xt=x; xt(i)=xt(i)+2.0d0*eps; xt(j)=xt(j)+2.0d0*eps
            call this%pot(xt,fb)
            xt=x; xt(i)=xt(i)+eps; xt(j)=xt(j)+eps
            call this%pot(xt,fc)
            xt=x; xt(i)=xt(i)-eps; xt(j)=xt(j)-eps
            call this%pot(xt,fd)
            xt=x; xt(i)=xt(i)-2.0d0*eps; xt(j)=xt(j)-2.0d0*eps
            call this%pot(xt,fe)
            hess(i,j)=(-fb+16.0d0*fc-30.0d0*fa+16.0d0*fd-fe)/(24.0d0*eps**2) - (hess(i,i)+hess(j,j))/2.0d0
            hess(j,i) = hess(i,j)
          end do
        end do
end subroutine

subroutine pes_hess(this, x, eps, hess)    !calculate hessian
        ! M. Herbst:
        !NOTE: We assume the correct coordinate convention is used and input is
        !in atomic units!
        class(pes), intent(in) :: this                                      ! selected pes
        double precision,dimension(9*this%Nw),intent(in) :: x               ! coordinates
        double precision, intent(in) :: eps                                 ! epsilon for hessian calc
        double precision,dimension(9*this%Nw,9*this%Nw),intent(out) :: hess ! hessian
        
        SELECT CASE(this%pes_type)
                case(PES_TTM3_F)
                        call pes_standard_hessian_from_potential2(this, x, eps, hess)
                case(PES_MBPOL)
                        call pes_five_point_stencil_hessian_from_potential(this, x, eps, hess)
                !ADD HERE
        end select
end subroutine

subroutine pes_pes_name(this,pes_string)       ! get the name of pot object
        class(pes), intent(in) :: this                   !selected pes
        character(len=8), intent(out) :: pes_string
        call pes_string_from_type(this%pes_type,pes_string)
end subroutine

subroutine pes_type_from_string(pes_string,pes_type)    ! get the pot object integer from string
        character(len=8), intent(in) :: pes_string
        integer,intent(out) :: pes_type

        select case (pes_string)
                case("TTM3-F")
                        pes_type = PES_TTM3_F
                case("MB-POL")
                        pes_type = PES_MBPOL
                !ADD HERE
                case default
                        pes_type = PES_NONE
        end select
end subroutine

subroutine pes_string_from_type(pes_type,pes_string)    ! get the pot name from pot object integer
        character(len=8), intent(out) :: pes_string
        integer,intent(in) :: pes_type

        SELECT CASE(pes_type)
                case(PES_TTM3_F)
                        pes_string = "TTM3-F"
                case(PES_MBPOL)
                        pes_string = "MB-POL"
                !ADD HERE
                case default
                        pes_string = "NONE"
        end select
end subroutine
end module
