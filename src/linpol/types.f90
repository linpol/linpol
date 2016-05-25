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
MODULE types
implicit none
!
! =======================================================================
! PURPOSE OF THIS MODULE
! =======================================================================
! Contains types to allow easy switching between single and double
! precision.
! 
!
! Record of revisions:
!
!  Date           Person         Change implemented
!  ----           ------         ------------------
!  ?? ??? ????    M. Cvitas      Development of module
!  22 Nov 2011    A. Reid        New version of module (based on MTC's)
!
! =======================================================================

      INTEGER,PARAMETER :: dp=kind(0.d0)
      INTEGER,PARAMETER :: sp=kind(0.0)

END MODULE types

