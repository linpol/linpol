#!/bin/bash
# This file is part of LINPOL, Copyright 2013-2016 Michael F. Herbst
# (info@michael-herbst.com)
#
# LINPOL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License can be found 
# in the file "LICENCE". If this file is not present in
# this package, see <http://www.gnu.org/licenses/>.
#


V="$2"		#= -v implies verbose
DIR="$1"

if [ "$1" == "-h" ]; then
	cat << EOF
	`basename $0` [ -h ] | <DIR> [ -v ]

Checks if the calculation in the directory <DIR> was successful.
Returns 0 if yes and 1 otherwise. If you want to see more details
provide a -v as second argument.
EOF
fi

if [ ! -d "$DIR" ]; then
	echo "Please provide a valid directory for the check, not $DIR"
	exit 1
fi
cd "$DIR"

#--------------------------------------
LOG=`ls -1t ????????-??????-log.txt | head -n1`
#           20121103-224530
[ "$V" == "-v" ] && echo -n "Is logfile present?  "
[ -z "$LOG" ] && { [ "$V" == "-v" ] && echo "...  no"; exit 1; }
[ "$V" == "-v" ] && echo "...  yes"
#--------------------------------------

#--------------------------------------
[ "$V" == "-v" ] && echo -n "Is L-BFGS converged?  "
< $LOG grep -s "L-BFGS CONVERGED SUCCESSFULLY" > /dev/null
[ "$?" != "0" ] && { [ "$V" == "-v" ] && echo "...  no"; exit 1; }
[ "$V" == "-v" ] && echo "...  yes"
#--------------------------------------

#--------------------------------------
[ "$V" == "-v" ] && echo -n "Is linpolexe completed?  "
< $LOG grep -s "Code complete" > /dev/null
[ "$?" != "0" ] && { [ "$V" == "-v" ] && echo "...  no"; exit 1; }
[ "$V" == "-v" ] && echo "...  yes"
#--------------------------------------

#--------------------------------------
[ "$V" == "-v" ] && echo -n "Is potential sensible?  "
RES=`< $LOG grep -A3 "PROPERTY CALCUALTION" | grep -o "potential:[[:space:]]*[^[:space:]]*" | awk '{if ($2 < 0) a=-$2; else a=$2; print (a < 1e-3)}'`
[ "$RES" == "1" ] && { [ "$V" == "-v" ] && echo "...  no:  too small"; exit 1; }
[ "$V" == "-v" ] && echo "...  yes"
#--------------------------------------

#TODO: check that a calculation has been finished correctly and no great errors have occurred

exit 0
