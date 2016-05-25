#!/bin/bash
# This file is part of LINPOL, Copyright 2013-2016 Michael F. Herbst
# (info@michael-herbst.com), Adam Reid & Jeremy Richardson
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

. linpol_scripts_common.sh

if [ "$1" == "-h" ]; then
	cat << EOF
	`basename $0` [ -h | <XYZ-File> <DoubleXYZ File> ]
Doubles the number of beads keeping beta the same.
EOF
	exit 0
fi

if [ ! -f "$1" ]; then
	echo "Please provide an Input xyz file (which contains beads to be doubled)"
	exit 1
fi

if [ -z "$2" ]; then
	echo "Please provide an Output xyz file (which contains the beads after doubling)"
	exit 1
fi

function doubling() {
	python << EOP
	#!/usr/bin/python
	#
	#Idea: Normal modes contain structural information of RING polymer (like where CofM is, the size in various directions, ...)
	#Hence if we double the beads in a way such that we keep the structure in normal mode coords as closely as possible
	#We preserve as much structure as possible.
	#Hence:

	import math, numpy, sys
	import vmd, nm

	NM = nm.cyclic()

	x, atomlist = vmd.load("$1", getatomlist=True)
	x = numpy.concatenate((x,x[::-1]))		#The method only works with cyclic polymer => make cyclic one by going forwards & backwards, ie pasting the linpol forwards and backwards together
	N = len(x)

	q = NM.forward(x)				#forward FT to get normal mode coords
	q = numpy.concatenate((q,numpy.zeros_like(q)))	#double coordinates by adding zeros for new normal mode coords (hence preserve strucutre of old)
	q *= math.sqrt(2)				#we doubled the coords => normalisation needs to be adapted by multiplication by sqrt(2)
	x = NM.backward(q)				#go back to real space

	beta = str("$3")

	vmd.xyz(x[:N], atomlist, "$2", beta)	#only keep first N beads (hence go back from ring back to linear polymer)
	EOP
}

#determine beta:
BETA=`cat "$1" | awk 'NR == 2 {print $1}'`

doubling "$1" "$2" "$BETA"
exit $?









