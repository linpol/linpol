#!/bin/bash
# This file is part of LINPOL, Copyright 2013-2016 Michael F. Herbst
# (info@michael-herbst.com).
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
	`basename $0` [ -h | <XYZ-File> <INDEX1> <INDEX2> ]
Permutes the two atoms at the 1based indices <INDEX1> and <INDEX2>.
Note that oxygens and hydrogens are refered to by numbers and this 
script does no checks whatsoever if you are not doing something
silly (ie swap hydrogen and oxygen or so)
EOF
	exit 0
fi

if [ ! -f "$1" ]; then
	echo "Please provide an existing xyz file containing the beads, not $1"
	exit 1
fi

if [ -z "$2" -o -z "$3" ]; then
	echo "Please provide two indices of the atoms to swap."
	exit 1
fi

if [ "$2" -eq "$3" ]; then
	echo "The indices should be different!"
	exit 1
fi

if [ "$2" -lt "$3" ]; then
	I1="$2"
	I2="$3"
else
	I1="$3"
	I2="$2"
fi

#Determine the number of atoms per bead
NA=`head -n1 "$1"`
#We need to change lines 2, (NA+2)+1, 2*(NA+2)+1

TMP=`mktemp`
cat "$1" | awk -v "na=$NA" -v "i1=$I1" -v "i2=$I2" 'NR == c*(na+2)+2+i1 { buffer=$0; next }; {print $0 }; NR == c*(na+2)+2+i2 {print buffer;c++}' > $TMP

mv $TMP "$1"
exit 0
