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
MODULE parameters
use clusterdatamod
use pes_interface
use phys_const          !the standard physical constants
implicit none
!
! For physical constants see phys_const.f90
! For parameters needed for the min module, see parameters_min.f90
!
!========================================!
!  Parameters read from command line     !
!========================================!
CHARACTER(len=64) :: filename                           ! XYZ File to open (possibly contains some beta information)
double precision :: masterbeta = 0.0                       ! If > 0 gets preference over the beta read in the input file (see prep.f90, xyzread for details)
INTEGER :: lbfgs_maxiter = 10000
logical :: debug = .false.
INTEGER :: endchoice = 2                                ! Choice of linear polymer: (1) fixed ends, (2) unfixed ends
double precision :: starteps = 1.0D-7                      ! Convergence criterion for L-BFGS search for single cluster geometry
double precision :: choseneps = 5.0D-5                     ! Convergence criterion for L-BFGS search for linear polymer geometry
double precision :: gradeps = 1.0E-3                         ! for Bowman gradient, whis is currently calculated as (f(x-eps)+f(x+eps))/2eps
                                                        ! hence numerical problems might occurr for smaller values
                                                        ! In a small number of tests best convergence was established for gradeps between 1E-3 and 1E-4
double precision :: hesseps = 2.0D-2                         ! for Bowman hessian, which is currently calculated via fhess2 (five-point stencil)
integer :: pes_type = PES_NONE                          ! which potential shall be used for the calculations 
type(pes) :: pes_obj                                    ! object representing the pes to use 
!
! for the selected opt_type see opt_selected.f90

!Controls which part of the program are run:
logical :: do_iter = .true.        ! should lbfgs iteration and optimisation be done
logical :: do_hess = .true.        ! should hessian be calculated and diagonalised
logical :: do_split = .true.       ! should Splitting parameters (logdet, phi, h) be calculated

!Reading and writingdata:
character(len=512) :: datainput = ''       ! where is the cluster data read from
character(len=512) :: dataoutputfile = ''  ! if single bead mode and nonempty, converged data is written to this file

!=======================================!
!  Other tweaks - recompile for effect  !
!=======================================!
integer :: writeline = 1        ! Only write L-BFGS output to screen for iterations that are multiples of this number
integer :: savestate = 10       ! Only save L-BFGS geometry to file for iterations that are multiples of this number

!==============================!
!  Internal runtime variables  !
!==============================!
INTEGER :: root = 0                                      ! Identity of root process number for MPI
!INTEGER :: mpi_n = 1                                    ! Number of mpi processes                                            TODO: implement
!INTEGER :: mpi_myrank = 1                               ! Rank of this mpi process                                           TODO: implement

!read in clusterdata:
type(clusterdata) :: cdata

!latest status in lbfgs:
double precision :: latest_gnorm                        ! contains the latest gnorm of lbfgs

!output files:
CHARACTER(len=64) :: initoutputfile = "unoptimized.xyz" ! File to create for unoptimized linear polymer
CHARACTER(len=64) :: outputfile = "optimized.xyz"       ! File to create for optimized linear polymer
CHARACTER(len=64) :: logfile = "log.txt"                ! File to log output
CHARACTER(len=64) :: pathwayfile = "pathway.txt"        ! File to save instanton pathway
CHARACTER(len=64) :: savefile = 'saved.xyz'             ! File to save partially-optimized linear polymer geometry
CHARACTER(len=64) :: lowestgradfile = "lowestgrad.xyz"  ! File to save linear polymer geometry with lowest gradient
CHARACTER(len=64) :: lowestfuncfile = "lowestfunc.xyz"  ! File to save linear polymer geometry with lowest function value











logical :: singleBeadMode = .false.!TODO: Temporary

END MODULE parameters
