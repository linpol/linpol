#!/bin/bash
#########################################
# LINPOL: LINear POLymer instanton code #
#########################################
#
# Research conducted using LINPOL or any derivative thereof should cite the following PhD thesis:
#
# A. A. Reid, "Quantum Tunnelling Splittings in Water Clusters, from Ring-Polymer Instanton Theory", University of Cambridge (2014)
#
#########################################
#
# This file is part of LINPOL, Copyright (C) Adam Reid, Jeremy Richardson &
# Michael F. Herbst 2011-2016.
#
# LINPOL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LINPOL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with LINPOL.  If not, see <http://www.gnu.org/licenses/>.
#
#########################################

#assumptions:	- mail has been installed
#		- linpolexe is available from the PATH variable		(alternatively set PAYLOAD according to your needs)
#
#By default mail is sent to the local inbox. 
#This can be changed by creating a file "~/.forward" and placing your preferred email address in there.

#------------------------------------------
#Settings:
PAYLOAD=""
FILENAME="complete.txt"
#------------------------------------------
#Runtime vars:
TO=$USER
PROC=1
[ -f /usr/bin/lscpu ] && PROC=`lscpu | grep "^CPU(s):" | awk '{print $2}'`
ARGS=

#------------------------------------------
prHelp() {
	cat << EOF
	`basename $0` [ -h | -x <Payload> | -t <TO@example.com> | -p <No Procs> | -- <ARGS for PAYLOAD> ]
Runs the payload ($PAYLOAD) using mpirun and emails the user once its done.

Options:
	-h		Displays this help.
	-t <TO>		Sends the email to the address <TO> instead of the current user
			Note that "null" disables sending mail
	--nomail	same as "-t null"
	-p <No Procs>	Uses this number of processes (Default: Number of processors on 
			current machine, or 1 if script can't determine that. 
			For your system the default would be $PROC)
	-x <Payload>	The payload program
	--		all following strings and arguments are passed straight onto the 
			payload program.
	other		All unrecognised switches and parameters are also passed onto
       			the payload program
EOF
}

#------------------------------------------

while [ "$1" ]; do
	case "$1" in
		"-h")	prHelp
			exit 0
			;;
		"-t")	shift
			TO="$1"
			;;
		"--nomail")
			TO="null"
			;;
		"-p")	shift
			PROC="$1"
			;;
		"-x")	shift
			PAYLOAD="$1"
			;;
		"--")	shift
			while [ "$1" ]; do
				ARGS="$ARGS \"$1\""
				shift
			done
			break
			;;
		*)	ARGS="$ARGS \"$1\""
			#TODO: !!!!!!!!!!
			;;
	esac
	shift
done

#-----------------------------------------

if [ -z "$PAYLOAD" ]; then 
	echo "No payload given"
	exit 1
fi

echo "-----------------------------------------------------------------------"
echo "Email will be sent to $TO upon completion"
echo "-----------------------------------------------------------------------"
echo

D1=`date`
#let's go:
eval "mpirun -np $PROC $PAYLOAD $ARGS"
RES=$?
D2=`date`

EXTRA="successfully"
[ "$RES" != "0" ] && EXTRA="with an ERROR !!"

( cat > $FILENAME ) << EOF
The program $PAYLOAD has finished $EXTRA!

Return code:	$RES
Host:		$HOSTNAME
PWD:		$PWD
Args:		$ARGS
Start:		$D1
End:		$D2

Hope all is well. ;)
EOF

if [ "$TO" == "null" ];then
	echo
	echo "-----------------------------------------------------------------------"
	echo "Done. Mail is not sent."
	echo "-----------------------------------------------------------------------"
	exit $RES
fi

echo
echo "-----------------------------------------------------------------------"
if type mail >> /dev/null; then
	#ie we have the mail program on the system
	cat $FILENAME | mail -s "[mpigo.sh]  $PAYLOAD finished $EXTRA" $TO

	echo "Done. Mail sent to $TO"
	RES=0
else
	echo "Done. Mail, however could not be sent !!!"
	RES=1
fi

echo "-----------------------------------------------------------------------"
exit $RES
