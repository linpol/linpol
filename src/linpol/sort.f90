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
!# This file is part of LINPOL, Copyright (C) Adam Reid, Jeremy Richardson & Michael Herbst 2011-2013.
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
module sort
implicit none
private         !everything private

public :: dmerge,dmergesort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- Merge Sort --!
!!!!!!!!!!!!!!!!!!

! based on the algo given at http://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran

contains
subroutine dmerge(n1,from1,n2,from2,nto,to)
        !merges two allocated, sorted arrays from1,from2 into an allocated array to
        double precision, dimension(n1),intent(in) :: from1
        double precision, dimension(n2),intent(in) :: from2
        double precision, dimension(nto),intent(inout) :: to
        integer,intent(in) :: n1, n2, nto
        integer :: i1,i2,ito
        !
        i1 = 1
        i2 = 1          !integers
        ito = 1
        !
        do while ( i1 <= n1 .and. i2 <= n2)
                if (from1(i1) <= from2(i2)) then
                        to(ito) = from1(i1)
                        i1 = i1+1
                else
                        to(ito) = from2(i2)
                        i2 = i2+1
                endif
                ito = ito+1
                !
                if ( ito > nto ) return ! to array is full
        end do
        !
        ! either with from1 or from2 we reached the end AND to array is not yet
        do while (i1 <= n1) !if for from1 we reached the end we will simply jump over this
                to(ito) = from1(i1)
                i1 = i1+1
                ito = ito+1
                if ( ito > nto ) return
        end do
        !
        do while (i2 <= n2)
                to(ito) = from2(i2)
                i2 = i2+1
                ito = ito+1
                if ( ito > nto ) return
        end do
        
        if ( ito <= nto) stop "Error in dmerge: size(from1)+size(from2) < size(to)"
end subroutine
! 
recursive subroutine dMergeSort(nA,A,T)     !only for internal usage; has no support for allocatable arrays
        !sorts the allocated array A, containing doubles
        !if T is present it is assumed to be allocated to (size(A)+1)/2 elements
        double precision, dimension(nA), intent(in out) :: A
        double precision, dimension( (nA+1)/2 ), intent(out) :: T
        integer, intent(in) :: nA
        double precision :: buffer

        integer :: i 
        ! Note: i is the index of division, ie all elements upt i go
        ! into first half, all following ones go into the second half

        if (nA < 2) return !nothing to do
        if (nA == 2) then
                if (A(1) > A(2)) then
                        !just swap them
                        buffer = A(1)
                        A(1) = A(2)
                        A(2) = buffer
                endif
                return !done for this case
        endif

        i=(nA+1)/2
        call dMergeSort(i,a,T)
        call dMergeSort(nA-i,a(i+1),T)

        if (A(i) > a(i+1)) then
                !merge them
                T(1:i)=A(1:i)
                call dmerge(i,T,nA-i,A(i+1),nA,A)
        endif
        return
end subroutine dMergeSort
end module
