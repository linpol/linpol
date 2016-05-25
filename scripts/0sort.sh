#!/bin/bash
# This file is part of LINPOL, Copyright 2013-2016 Michael F. Herbst
# (info@michael-herbst.com), Jeremy Richardson & Adam Reid
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

. linpol_scripts_common.sh

if [ "$1" == "-h" ]; then
	cat << EOF
	`basename $0` [ -h | <XYZ-File> ]
Genereates a <XYZ-File>_sorted.xyz and a 
<XYZ-File>_labelled.xyz containing the output.
EOF
	exit 0
fi

if [ -z "$1" ]; then
	echo "Please provide an xyz file"
	exit 1
fi

function sorting() {
	# the python bit
python << EOP

# =======================================================================
# PURPOSE OF THIS PROGRAM
# =======================================================================
# This program was written to analyse a water cluster geometry in order
# to match each oxygen atom to a pair of hydrogen atoms.  It outputs two
# xyz files, each with the lines ordered to match the requirements of the
# Bowman potential code:
#
# 'Atoms need to be in the order of ascending mass, so "H" always goes
# before "O". All "H" atoms and "O" atoms should be listed in the same
# monomer order.  For example, a water dimer should look like, "H1 H2 H3
# H4 O5 O6", where H1-O5-H2 is one monomer and H3-O6-H4 is the second
# monomer.'
# (Bowman readme)
#
# The 'labelled' output file labels each trio of atoms with the same
# letter of the alphabet (so the analysis can be checked visually in VMD).
#
# The 'unlabelled' output file simply uses the usual 'H' and 'O' labels,
# and hence is appropriate for use as a starting geometry in the linpol
# code.
#
# Command line arguments:
# The first argument is the input file, the second is the name for the
# 'labelled' output file and the third is the name for the 'unlabelled'
# output file.
#
# Record of revisions:
#
#  Date           Person         Change implemented
#  ----           ------         ------------------
#  4 Jul 2012     A. Reid        Original code developed
#
# =======================================================================


########################################################################################################################
# Import modules and functions
########################################################################################################################

import sys                        
from linpolmodule import datafile
import string
from math import sqrt
from numpy import zeros
from array import array

########################################################################################################################
# Set initial parameters and files
########################################################################################################################

inputfile = "$1"
outputfilelabelled = "$2"
outputfileunlabelled = "$3"

########################################################################################################################
# Initial messages
########################################################################################################################

print('==========================================================================')
print('WATER GEOMETRY ANALYSIS CODE')
print('==========================================================================')
print('')
print('***Input file selected by user***') 
print('Input file: ' + str(inputfile))
print('')

########################################################################################################################
# Set output file
########################################################################################################################

outputlabelled = open(str(outputfilelabelled) , 'w')
outputunlabelled = open(str(outputfileunlabelled) , 'w')

########################################################################################################################
# Print details of input file
########################################################################################################################

geometry = datafile(inputfile)

print('***Analyzing input file***')
nlines = geometry.readlines()
print('Number of lines in file: ' + str(nlines))
natoms = geometry.readatoms()
print('Number of atoms in file: ' + str(natoms))
nwaters = natoms/3
print('Number of water molecules in file: ' + str(nwaters))
beta = geometry.readbeta()
print('Beta in file: ' + str(beta))
print('')
if nwaters > 26:
    print('Atom labels dictionary has only 26 entries.')
    print('To use more than 26 water molecules, extend the dictionary.')
    print('')
    sys.exit('***SCRIPT ENDING PREMATURELY: Dictionary failure***')

########################################################################################################################
# Set up labels dictionary and distance function
########################################################################################################################

labelsdict = {}
index = 1
for letter in string.ascii_uppercase:
    labelsdict[index] = letter
    index += 1

def calcdist( Oindex, Hindex ):
    xdiff = x[1,Oindex] - x[1,Hindex]
    ydiff = x[2,Oindex] - x[2,Hindex]
    zdiff = x[3,Oindex] - x[3,Hindex]
    dist = sqrt(xdiff**2 + ydiff**2 + zdiff**2)
    return dist

########################################################################################################################
# Analyse geometry
########################################################################################################################

print('***Analyzing water cluster geometry***')
x = geometry.readall()
x= geometry.x
xout = zeros((4,natoms))

distances = zeros(natoms)	# Blank array for distances

molecule = 1
entrycount = 0

for entry in x[0]:

    if entry == 1:
	Oindex = entrycount

	other = 0
	while other < natoms:
	    distances[other] = calcdist(Oindex,other)
	    other += 1

	order = distances.argsort()
	Honeindex = order[1]
	Htwoindex = order[2]

	outOindex = ((natoms*2)/3) + (molecule-1)
	outHtwoindex = (molecule*2)-1
	outHoneindex = outHtwoindex-1

	xout[0,outOindex] = molecule
	xout[0,outHoneindex] = molecule
	xout[0,outHtwoindex] = molecule

	for coord in [1,2,3]:
	    xout[coord,outOindex] = x[coord,Oindex]
	    xout[coord,outHoneindex] = x[coord,Honeindex]
	    xout[coord,outHtwoindex] = x[coord,Htwoindex]

	molecule += 1

    entrycount += 1

########################################################################################################################
# Write to labelled output file
########################################################################################################################

outputlabelled.write(str(natoms) + '\\n')      # No. of atoms
outputlabelled.write(str(beta) + '\\n')        # Beta
    
atom = 0    
while atom < natoms:
    entry = str(labelsdict[xout[0,atom]]) + ' ' + str(xout[1,atom]) + ' ' + str(xout[2,atom]) + ' ' + str(xout[3,atom]) + '\\n'
    outputlabelled.write(entry.replace('e','E'))
    atom = atom + 1

print('The labelled geometry has been written to: ' + str(outputfilelabelled)) 

########################################################################################################################
# Write to unlabelled output file
########################################################################################################################

outputunlabelled.write(str(natoms) + '\\n')      # No. of atoms
outputunlabelled.write(str(beta) + '\\n')        # Beta
    
atom = 0
while atom < natoms:
    if atom >= (natoms*2)/3:
	label = 'O'
    else:
	label = 'H'
    entry = str(label) + ' ' + str(xout[1,atom]) + ' ' + str(xout[2,atom]) + ' ' + str(xout[3,atom]) + '\\n'
    outputunlabelled.write(entry.replace('e','E'))
    atom = atom + 1

print('The unlabelled geometry has been written to: ' + str(outputfileunlabelled)) 
print('')

########################################################################################################################
# Close files
########################################################################################################################

outputlabelled.close()
outputunlabelled.close()

########################################################################################################################
# Closing messages
########################################################################################################################

print('==========================================================================')
print('Done')
print('==========================================================================')

EOP
}

INPUT="$1"
OUTPUT="`basename $1 .xyz`_sorted.xyz"
LABELLED="`basename $1 .xyz`_labelled.xyz"
TMP=`mktemp`

#Adams Skript requires a beta which I think should not be neccessary => add one in case it's not there and remove later
LINE2=`cat $1 | awk 'NR == 2 {print; exit}'`
cat $INPUT | awk 'NR == 2 {print "0"; next}; {print $0}' > $TMP
sorting "$TMP" "$LABELLED.tmp" "$OUTPUT.tmp"
RES=$?

< "$LABELLED.tmp" awk -v "line=$LINE2" 'NR == 2 {print line; next}; {print $0}' > $LABELLED && rm $LABELLED.tmp
< "$OUTPUT.tmp" awk -v "line=$LINE2" 'NR == 2 {print line; next}; {print $0}' > $OUTPUT && rm $OUTPUT.tmp

rm -r $TMP
exit $RES
