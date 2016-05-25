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

if [ "$1" == "-h" ]; then
	cat << EOF
	`basename $0` [ -h | <XYZ-File> <BETA> ]
Changes the beta value in the <XYZ-File> to BETA for all beads.
EOF
	exit 0
fi

if [ ! -f "$1" ]; then
	echo "Please provide an existing xyz file containing the beads, not $1"
	exit 1
fi

if [ -z "$2" ]; then
	echo "Please provide a new beta value."
	exit 1
fi

#Determine the number of atoms per bead
NA=`head -n1 "$1"`
#We need to change lines 2, (NA+2)+1, 2*(NA+2)+1

TMP=`mktemp`
cat "$1" | awk -v "na=$NA" -v "beta=$2" 'NR == c*(na+2)+2 { print beta; c++; next }; {print $0 }' > $TMP

mv $TMP "$1"
exit 0
