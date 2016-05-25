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
module clusterdatamod
use pes_interface
implicit none
private

public :: is_valid_cluster_data_string, get_cluster_data_full_path, readClusterData, clusterdata, writeClusterData

type :: clusterdata
      integer :: dof = 0 !degrees of freedom = d*p
      double precision :: gnorm = 0d0
      character(len=15) :: dateTime = ''
      logical :: global_minimum = .false.
      character(len=512) :: comments = ''
      double precision :: potential = 0d0
      double precision, allocatable,dimension(:) :: evalues
      character(len=512) :: fullPath = '' !not written if --writedata; for internal usage
end type

interface writeClusterData
        module procedure writeClusterDataContainer !writes clusterdata type to a file
        module procedure writeClusterDataCalled ! writes some things passed to the function to a file
end interface

contains
function get_cluster_data_prefix(pes_obj)
        character(len=512) :: get_cluster_data_prefix
        CHARACTER(LEN=512) :: dir
        CHARACTER(LEN=8) :: pesN
        type(pes),intent(in) :: pes_obj
        CALL GET_ENVIRONMENT_VARIABLE("LINPOL", dir)
        call pes_obj%pes_name(pesN)
        get_cluster_data_prefix = trim(dir)//'/data/'//trim(pesN)//'/cluster'
end function

function get_cluster_data_default()
        character(len=512) :: get_cluster_data_default
        CHARACTER(LEN=512) :: tmp
        CALL GET_ENVIRONMENT_VARIABLE("LINPOL_DEFAULT_CLUSTER", tmp)
        get_cluster_data_default = trim(tmp)
end function

function get_cluster_data_full_path(instr,pes_obj,success)
        !returns full path and success = true if it works, otherwise false and
        !empty string

        logical,intent(out) :: success !did it work at all
        character(len=512),intent(in) :: instr ! the string to expand to the full path
        type(pes),intent(in) :: pes_obj
        character(len=512) :: str
        character(len=512) :: get_cluster_data_full_path

        integer :: str_len
        CHARACTER(80) :: msg                    ! Status message
        INTEGER :: stat                         ! IO status message

        ! if str does not end with .dat, add it:
        str_len = len_trim(instr)
        if (str_len .gt. 3) then
                if ( instr(str_len-3:str_len) .ne. '.dat' ) then
                        str=trim(instr)//".dat"
                else
                        str=instr
                end if
        else
                !string is too short to incorporate the .dat ending!
                str=trim(instr)//".dat"
        end if

        success = .true.
        !===========================================================================

        ! 1. Try to simply open what's given:
        OPEN (UNIT=999, FILE=str, STATUS='OLD', ACTION='READ', IOSTAT=stat, IOMSG=msg )
        if (stat .eq. 0) then
                get_cluster_data_full_path = str
                return
        end if
        CLOSE (UNIT=999)

        ! 2. Try to add data prefix at the front
        get_cluster_data_full_path = trim(get_cluster_data_prefix(pes_obj))&
                //'/'//trim(str)
        OPEN (UNIT=999, FILE=get_cluster_data_full_path, STATUS='OLD', ACTION='READ', IOSTAT=stat, IOMSG=msg )
        if (stat .eq. 0) return
        CLOSE (UNIT=999)

        ! 3. Try to add default in the middle & prefix at front
        get_cluster_data_full_path = trim(get_cluster_data_prefix(pes_obj)) &
                //'/' // TRIM(get_cluster_data_default()) //'/'// trim(str)
        OPEN (UNIT=999, FILE=get_cluster_data_full_path, STATUS='OLD', ACTION='READ', IOSTAT=stat, IOMSG=msg )
        if (stat .eq. 0) return
        CLOSE (UNIT=999)

        ! All went wrong:
        get_cluster_data_full_path = ''
        success = .false.
end function

function is_valid_cluster_data_string(str,pes_obj) result(valid)
        logical :: valid
        character(len=512) :: str
        type(pes) :: pes_obj
        character(len=1) :: outp
        outp = get_cluster_data_full_path(str,pes_obj,valid)
end function

subroutine readClusterData(this, str,pes_obj, success)
        ! try to read cluster data. Success is true if it works and this is filled
        ! with the data or false if it went wrong and this is empty.

        type(clusterdata),intent(out) :: this
        character(len=512),intent(in) :: str
        type(pes),intent(in) :: pes_obj
        logical, intent(out) :: success
        character(len=512) :: path

        character(len=512) :: line
        INTEGER :: stat
        logical :: potential_read = .false. !has the potential been read (required for correct run of program!)

        path = get_cluster_data_full_path(str,pes_obj,success)
        if (.not. success) then
                return
        end if
        
        ! open file
        OPEN (UNIT=999, FILE=path, STATUS='OLD', ACTION='READ', IOSTAT=stat )
        if (stat .ne. 0) then
                success = .false.
                return
        end if
        this%fullPath = path

        !read loop
        do while (.true.)
               989 continue
               READ ( 999, '(A)', IOSTAT=stat ) line 

               if (stat .lt. 0) then
                        !end of file reached => check if all required data has
                        !been read

                        if ( (this%dof .gt. 0 .and. potential_read) .and. allocated(this%evalues) ) exit
                                
                        !some things are missing
                        success = .false.
                        CLOSE (UNIT=999)
                        return
                end if

                if (stat .gt. 0 ) then
                        !some error occurred:
                        success = .false.
                        CLOSE (UNIT=999)
                        return
                end if

                ! --------------------------- COMMENT LINES
                if (line(1:1) .eq. '!' .or. line(1:1) .eq. '#' ) then
                        !comment line
                        goto 989 !go to start of loop
                end if

                !---------------------------- EVALUES

                if (trim(line) .eq. "-------------------------------START EVALUES-----------------------------") then
                        !read evalues:

                        if (this%dof .le. 0) then
                                ! this is a requirement and needs to be read
                                ! befor we read this
                                success = .false.
                                CLOSE (UNIT=999)
                                return
                        end if

                        allocate(this%evalues(this%dof))
                        read(999, *, IOSTAT=stat) this%evalues
                        if (stat .ne. 0) then
                                !some error occurred:
                                success = .false.
                                CLOSE (UNIT=999)
                                return
                        end if
                        
                        READ ( 999, '(A)', IOSTAT=stat ) line
                        
                        if (trim(line) .ne. "--------------------------------END EVALUES------------------------------") then
                                !not correct format
                                success = .false.
                                CLOSE (UNIT=999)
                                return
                        end if

                        !read next line:
                        goto 989 !start of loop
                end if

                ! ----------------- FIELDS
                ! line can only be a field                 

                if (line(21:23) .ne. ' = ') then
                        !that should not be the case
                        success = .false.
                        CLOSE (UNIT=999)
                        return
                end if

                !Required to not have any crap in the variable:
                this%comments=''
                select case(trim(line(1:20)))
                        case('dof')
                                READ(line(24:),*,IOSTAT=stat) this%dof
                                if (stat .ne. 0) then
                                        success = .false.
                                        CLOSE (UNIT=999)
                                        return
                                end if
                        case('gnorm')
                                READ(line(24:),*,IOSTAT=stat) this%gnorm
                                if (stat .ne. 0) then
                                        success = .false.
                                        CLOSE (UNIT=999)
                                        return
                                end if
                        case('date-time')
                                this%dateTime=TRIM(line(24:))
                        case('global_minimum')
                                READ(line(24:),*,IOSTAT=stat) this%global_minimum
                                if (stat .ne. 0) then
                                        success = .false.
                                        CLOSE (UNIT=999)
                                        return
                                end if
                        case('comments')
                                this%comments=trim(line(24:))
                        case('potential')
                                READ(line(24:),*,IOSTAT=stat) this%potential
                                if (stat .ne. 0) then
                                        success = .false.
                                        CLOSE (UNIT=999)
                                        return
                                end if
                                potential_read=.true.
                        case default
                                !error
                                success = .false.
                                CLOSE (UNIT=999)
                                return
                end select
        end do

        CLOSE (UNIT=999)
        success = .true.
        return
end subroutine

! Specification for the cluster dat files
! FIELD            = VALUE
! <-- 20 Chars --><3>
!
! NOTE: The field has to be 20 chars long, followed by " = " and the value of a
! format, which can be parsed to the required data type.
!
! Possible FIELDs
!       - dof (integer)         degrees of freedom
!       - gnorm (double)        the gnorm of the structure for this data
!       - date-time (char15)    date & time
!       - global_minimum (logical)   is this the global min?
!       - comments (char512)    comments string (will be read in)
!       - potential (double)    the potential obtained
!
! the evalues in ASCENDING order, between
! -------------------------------START EVALUES-----------------------------
! and
! --------------------------------END EVALUES------------------------------
!
! the order of the fields is arbitrary, but the evalues block SHOULD be at the
! end and it MUST be after dof (which is required to read it in)
!
! Not all fields need to be specified, but the evalues, the potential and
! dof have to be.
!
! extra comment lines are permitted and start with ! or # and are not read in

function writeClusterDataContainer(filename, this, writeEverything) result(succ)
        type(clusterdata),intent(in) :: this
        character(len=512),intent(in) :: filename       ! the file to write to
        logical,optional,intent(in) :: writeEverything  ! write everything from this or only the things the program can know
                        ! after beeing successfully executed, ie 
                        !    - dof
                        !    - eps
                        !    - dateTime
                        !    - potential
                        !    - evalues

        logical :: succ 
        if (present(writeEverything) .and. writeEverything) then
                succ = writeClusterDataCalled(filename, this%dof, this%gnorm, &
                        this%dateTime, this%potential, this%evalues)
        else
                succ = writeClusterDataCalled(filename, this%dof, this%gnorm, &
                        this%dateTime, this%potential, this%evalues, &
                        this%global_minimum, this%comments)
        end if
end function

function writeClusterDataCalled(filestr, dof, gnorm, dateTime, potential, evalues, global_minimum, comments) result(succ)
        !return true if all was fine, else false

        character(len=512),intent(in) :: filestr
        integer,intent(in) :: dof
        double precision,intent(in) :: gnorm
        character(len=15),intent(in) :: dateTime
        double precision, intent(in) :: potential
        double precision, intent(in), dimension(dof) :: evalues

        logical,optional,intent(in) :: global_minimum
        character(len=512),optional,intent(in) :: comments

        character(len=512) :: filename
        character(len=20) :: field
        character(len=3) :: equals = ' = '
        integer :: stat
        integer :: str_len ! length of input string filestr
        logical :: succ
        succ=.true.

        str_len = len_trim(filestr)
        if (str_len .gt. 3) then
                if ( filestr(str_len-3:str_len) .ne. '.dat' ) then
                        filename=trim(filestr)//".dat"
                else
                        filename=filestr
                end if
        else if (str_len .eq. 0) then
                succ= .false.
                return
        else
                !string is too short to incorporate the .dat ending!
                filename=trim(filestr)//".dat"
        end if

        !Try to open the file for writing:
        OPEN (UNIT=999, FILE=filename, STATUS='REPLACE', ACTION='WRITE',IOSTAT=stat)

        if (stat .ne. 0) then
                CLOSE (UNIT=999)
                succ = .false.
                return
        end if
                
        ! write all given fields:
        field="dof"
        write(999,'(A,A,I20)') field, equals, dof

        field="gnorm"
        write(999,'(A,A,E20.13)') field, equals, gnorm
                
        field="date-time"
        write(999,'(A,A,A)') field, equals, trim(dateTime)

        field="global_minimum"
        if (present(global_minimum)) then
                write(999,'(2a,L20)') field, equals, global_minimum
        else
                write(999,'(3a,L20)') "!",field, equals, .false.
        end if

        field="comments"
        if (present(comments)) then
                write(999,'(3A)') field, equals, trim(comments)
        else
                write(999,'(4A)') "!",field, equals, ""
        end if

        field="potential"
        write(999,'(A,A,E20.13)') field, equals, potential

        ! Write evalues:
        write(999,'(A)') "-------------------------------START EVALUES-----------------------------"
        write(999,*) evalues
        write(999,'(A)') "--------------------------------END EVALUES------------------------------"

        CLOSE (UNIT=999)
end function
end module
