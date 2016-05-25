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
	`basename $0` [ -h | <XYZ-File> <HalfXYZ File> ]
Skript to throw away half of the beads in the xyz File.

Since many of the beads at the beginning and the end of the linear polymer
closely resemble the structure of the minima they connect, more beads are 
disposed of at the ends of the linear polymer:
	From the   1st         quarter only   1 in 4   beads is kept.
	From the   2nd & 3rd   quarter only   3 in 4   beads are kept.
	From the   4th         quarter only   1 in 4   beads is kept.
EOF
	exit 0
fi

if [ ! -f "$1" ]; then
	echo "Please provide an Input xyz file (which contains beads to be halfed)"
	exit 1
fi

if [ -z "$2" ]; then
	echo "Please provide an Output xyz file (which contains the beads after halfing)"
	exit 1
fi


#Determine the number of atoms per bead
NA=`head -n1 "$1"`
#Determine the number of beads:
NB=`cat "$1" | awk '
	BEGIN { line1="";line2="";stage=0; c=0 }
	NR == 1 {line1=$0; next}
       	NR == 2 {line2=$0;c++;next}; 
	$0 == line1 {stage=1;next};
	$0 == line2 && stage == 1 {c++;stage=0;next}	#increment if previous line was line1
	stage == 1 { stage=0 }				#previous line was line1, but this is not line2 => reset
	END{print c}'`
#hence one block of beads is NA+2 lines

cat "$1" | awk -v "na=$NA" -v "nb=$NB" \
	'BEGIN {
			pr=0		#do we print this block
			blockline=0	#the line in the current bead block (0 means disabled for this bead)
			bead=0		#current bead

			nblocks=4	#Number of blocks we partition the beads in
			keep[0]=1
			keep[1]=3	#fraction x/nblocks of beads to keep; expressed as x
			keep[2]=3
			keep[3]=1

			perBlock=int(nb/nblocks)
		}
	
		function getBlock(bead) {
			#returns the block number a bead is in.
			t=int((bead-1)/perBlock)	#bead is 1-based index
			return (t+1)
		}

		function keepBead(bead) {
			#returns 1 to keep the bead ; 0 to discard

			if (keep[getBlock(bead)-1] > (bead % nblocks)) {	#keep is 0-based index; getBlock is 1-based
				return 1
			}
			return 0
		}

		blockline==0 && keepBead(++bead)==1 {blockline=1;pr=1 }		#we print the bead
		blockline==0 {blockline=1;pr=0}					#we wont print the bead (increment done above)
		pr==1 {print $0}						#print the line
		blockline >= na+2 { blockline=0;next }				#this bead is finished => go to next
		blockline > 0 { blockline++}					#increment bead block line counter
' > "$2"
exit $?
