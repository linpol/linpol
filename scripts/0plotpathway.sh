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

IN=
OUT=

if [[ "$1" = "-h" || "$2" = "-h" || "$3" = "-h" ]]; then
	cat << EOF
	`basename $0` <PATHWAYFILE> [ <OUT> ]

Uses gnuplot to generate png plot of the pathwayfile given as <PATHWAYFILE>. 
If <OUT> is not given, the png file will have the same name as the input, 
but the extension replaced by ".png"

EOF
	exit 0
fi

if [ -d "$1" ]; then
	IN=`ls $1/????????-??????-pathway.txt 2> /dev/zero | head -n1` 
	if [ -z "$IN" ]; then
		echo "Error, invalid directory $1: Does not contain pathway file."
		exit 1
	fi
elif [ "$1" ]; then
	if [ ! -f "$1" ]; then
		echo "Error: This file or directory does not exist: $1"
		exit 1
	fi
	IN="$1"
else 
	IN=`ls ????????-??????-pathway.txt 2> /dev/zero | head -n1` 
	if [ -z "$IN" ]; then
		echo "Error, could not fine a pathway file in current directory."
		exit 1
	fi
fi

if [ "$2" ]; then
	OUT="$2"
else
	OUT="`basename "$IN" .txt`.png"
fi

/usr/bin/gnuplot << EOF
set terminal png size 900,600

set title "Pathway potential plot"
set xlabel "Mass weighted distance r"
set ylabel "Potential V / Hartree"

set output "$OUT"
plot "$IN" using 1:2 with linespoints
EOF

exit $?
