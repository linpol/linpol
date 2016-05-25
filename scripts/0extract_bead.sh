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

N="$2"
INPUT="$1"
OUT="$3"

if [ "$1" == "-h" ]; then
	cat << EOF
	`basename $0` [ -h | <XYZ-File> <N> <OUT> ]
Reads <XYZ-file> and prints the <N>th bead to <OUT>. If missing or equal to "-"
stdout is used. 

<N> can be
   number	the <N>th bead is extracted
   N		the last bead is extracted
   N/2		the halfway bead is extracted
EOF
	exit 0
fi

if [ ! -f "$INPUT" ]; then
	echo "Please provide an existing xyz file containing the beads, not $1" >&2
	exit 1
fi

if [ -z "$N" ]; then
	echo "Please provide a value for N" >&2
	exit 1
fi


#Determine number of beads:
NBEAD=`cat "$INPUT" | awk '
	BEGIN { line1="";line2="";stage=0; c=0 }
	NR == 1 {line1=$0; next}
	NR == 2 {line2=$0;c++;next}; 
	$0 == line1 {stage=1;next};
	$0 == line2 && stage == 1 {c++;stage=0;next}	#increment if previous line was line1
	stage == 1 { stage=0 }				#previous line was line1, but this is not line2 => reset
	END{print c}'`

#Determine the number of atoms per bead
NA=`head -n1 "$1"`

if [ "$N" = "N" ]; then
	N=$NBEAD
elif [ "$N" = "N/2" ]; then
	N=$((NBEAD/2))
elif [ $N -gt $NBEAD ];then
	echo "The given value for N is larger than the number of beads available!" >&2
	exit 1
fi

# a bead has NA+2 lines in the xyz file

BEAD=`< "$INPUT" awk -v "from=$(( (N-1)*(NA+2)+1 ))" -v "to=$(( N*(NA+2) ))" 'BEGIN {pr=0}; NR == from {pr=1}; pr==1; NR == to {pr=0}'`

if [ "$OUT" ] ; then
	echo "$BEAD" > "$OUT"
else
	echo "$BEAD"
fi

exit 0
