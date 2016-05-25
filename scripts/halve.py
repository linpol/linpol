#!/usr/bin/python
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

# =======================================================================
# PURPOSE OF THIS PROGRAM
# =======================================================================
# 
#
# Record of revisions:
#
#  Date           Person         Change implemented
#  ----           ------         ------------------
#  9 May 2013     A. Reid        Developed from build.py
#
# =======================================================================


########################################################################################################################
# Import modules and functions
########################################################################################################################

import sys                        
from numpy import zeros
from linpolmodule import datafile

def atomlabel(inputlinpol,atom,bead):
    if inputlinpol[0,atom,bead] == 0:
        label = 'H'
    elif inputlinpol[0,atom,bead] == 1:
        label = 'O'
    else:
        print('Error in atomic labels!')
    return label

########################################################################################################################
# Set initial parameters and files
########################################################################################################################

inputfile = str(sys.argv[1])
outputfile = str(sys.argv[2])

########################################################################################################################
# Initial messages
########################################################################################################################

print('==========================================================================')
print('INSTANTON STARTING GEOMETRY CREATION CODE')
print('==========================================================================')
print('')
print('***Files selected by user***') 
print('File to be halved: ' + str(inputfile))
print('New file to be created: ' + str(outputfile))
print('')

########################################################################################################################
# Set data files
########################################################################################################################

output = open(str(outputfile) , 'w')

########################################################################################################################
# Analysis of input files
########################################################################################################################

inputfileclass = datafile(inputfile)

natoms = inputfileclass.readatoms()
beta = inputfileclass.readbeta()
nlines = inputfileclass.readlines()
nbeads = nlines/(natoms+2)

inputlinpol = inputfileclass.readlinpol()

print('***Input file information***')
inputfileclass.writelines()
inputfileclass.writeatoms()
print('Hence, number of beads: ' + str(nbeads))
inputfileclass.writebeta()

########################################################################################################################
# Generate new linear polymer geometry
########################################################################################################################

print('')
print('***Generating new linear polymer geometry***')

bead = 0

while bead < nbeads:

    if bead/2 == bead/2.0:

        output.write(str(natoms) + '\n')      # No. of atoms
        output.write(str(beta) + '\n')        # Beta

        atom = 0    
        while atom < natoms:
            entry = str(atomlabel(inputlinpol,atom,bead)) + ' ' + str(inputlinpol[1,atom,bead]) + ' ' + str(inputlinpol[2,atom,bead]) + ' ' + str(inputlinpol[3,atom,bead]) + '\n'
            output.write(entry.replace('e','E'))
            atom = atom + 1

    bead = bead + 1

print('')
print('The new ' + str(nbeads/2) + '-bead linear polymer has been written to: ' + str(outputfile)) 
print('')

########################################################################################################################
# Close files
########################################################################################################################

output.close()

########################################################################################################################
# Closing messages
########################################################################################################################

print('==========================================================================')
print('Done')
print('==========================================================================')
