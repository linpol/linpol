#!/usr/bin/python
# Filename: build.py
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
# This program was written to calculate a new linear-polymer geometry
# (described by an xyz file) to act as a starting point for L-BFGS
# minimization.  It takes three input files: (1) minimized single bead
# geometry, (2) guess of halfway geometry, (3) tunnelled single bead
# geometry (similar to the minimized geometry, but with 1-2 atoms
# rearranged.  The fourth argument is the number of beads (N) and the
# the fifth is the name for the output file.
#
# An N-bead geometry is created, with the following numbers of beads
# from each input file: (1) (N/2)-(N/16) (2) N/8 (3) (N/2)-(N/16).  For
# example, when N=16, we have (1) 7 (2) 2 (3) 7, and when N=32 we have
# (1) 14 (2) 4 (3) 14.
#
# Record of revisions:
#
#  Date           Person         Change implemented
#  ----           ------         ------------------
#  28 Jun 2012    A. Reid        Original code developed
#  4 Jul 2012     A. Reid        Commented-out lines deleted
#
# =======================================================================


########################################################################################################################
# Import modules and functions
########################################################################################################################

import sys                        
from numpy import zeros
from linpolmodule import datafile

########################################################################################################################
# Set initial parameters and files
########################################################################################################################

firstbeadfile = str(sys.argv[1])
halfwaybeadfile = str(sys.argv[2])
lastbeadfile = str(sys.argv[3])
nbeads = int(sys.argv[4])
outputfile = str(sys.argv[5])

########################################################################################################################
# Initial messages
########################################################################################################################

print('==========================================================================')
print('INSTANTON STARTING GEOMETRY CREATION CODE')
print('==========================================================================')
print('')
print('***Input files selected by user***') 
print('First bead: ' + str(firstbeadfile))
print('Halfway bead: ' + str(halfwaybeadfile))
print('Last bead: ' + str(lastbeadfile))
print('')

########################################################################################################################
# Set data files
########################################################################################################################

output = open(str(outputfile) , 'w')

########################################################################################################################
# Analysis of input files
########################################################################################################################

firstbeadclass = datafile(firstbeadfile)
halfwaybeadclass = datafile(halfwaybeadfile)
lastbeadclass = datafile(lastbeadfile)

print('***Comparing files***')
if (firstbeadclass.readlines() == halfwaybeadclass.readlines() == lastbeadclass.readlines()):
    nlines = firstbeadclass.readlines()
    print('Same number of lines in each file?.. OK (' + str(nlines) + ')')
else:
    print('Same number of lines in each file?.. No!')
    firstbeadclass.writelines()
    halfwaybeadclass.writelines()
    lastbeadclass.writelines()
    sys.exit('***SCRIPT ENDING PREMATURELY: Input file failure***')
if (firstbeadclass.readatoms() == halfwaybeadclass.readatoms() == lastbeadclass.readatoms()):
    natoms = firstbeadclass.readatoms()
    print('Same number of atoms in each file?.. OK (' + str(natoms) + ')')
else:
    print('Same number of atoms in each file?.. No!')
    firstbeadclass.writeatoms()
    halfwaybeadclass.writeatoms()
    lastbeadclass.writeatoms()
    sys.exit('***SCRIPT ENDING PREMATURELY: Input file failure***')
if (firstbeadclass.readbeta() == halfwaybeadclass.readbeta() == lastbeadclass.readbeta()):
    beta = firstbeadclass.readbeta()
    print('Same beta in each file?.. OK (' + str(beta) + ')')
else:
    print('Same beta in each file?.. No!')
    firstbeadclass.writebeta()
    halfwaybeadclass.writebeta()
    lastbeadclass.writebeta()
    sys.exit('***SCRIPT ENDING PREMATURELY: Input file failure***')
if (firstbeadclass.readatomtags() == halfwaybeadclass.readatomtags() == lastbeadclass.readatomtags()):
    atomtags = firstbeadclass.readatomtags()
    print('Same elements in each file?.. OK')
else:
    print('Same elements in each file?.. No!')
    sys.exit('***SCRIPT ENDING PREMATURELY: Input file failure***')
print('')

########################################################################################################################
# Generate new linear polymer geometry
########################################################################################################################

print('***Generating new linear polymer geometry***')

nfirst = (nbeads/2)-(nbeads/16)
nhalfway = (nbeads/8)
nlast = (nbeads/2)-(nbeads/16)
print('Copies of first bead: ' + str(nfirst))
print('Copies of halfway bead: ' + str(nhalfway))
print('Copies of last bead: ' + str(nlast))

bead = 0

while bead < nbeads:

    output.write(str(natoms) + '\n')      # No. of atoms
    output.write(str(beta) + '\n')        # Beta
    
    if bead < nfirst:
        atom = 0    
        while atom < natoms:
            entry = str(firstbeadclass.readatomtags()[atom]) + ' ' + str(firstbeadclass.readcoords()[0,atom]) + ' ' + str(firstbeadclass.readcoords()[1,atom]) + ' ' + str(firstbeadclass.readcoords()[2,atom]) + '\n'
            output.write(entry.replace('e','E'))
            atom = atom + 1

    if nfirst <= bead < (nfirst + nhalfway):
        atom = 0    
        while atom < natoms:
            entry = str(halfwaybeadclass.readatomtags()[atom]) + ' ' + str(halfwaybeadclass.readcoords()[0,atom]) + ' ' + str(halfwaybeadclass.readcoords()[1,atom]) + ' ' + str(halfwaybeadclass.readcoords()[2,atom]) + '\n'
            output.write(entry.replace('e','E'))
            atom = atom + 1

    if (nfirst + nhalfway) <= bead:
        atom = 0    
        while atom < natoms:
            entry = str(lastbeadclass.readatomtags()[atom]) + ' ' + str(lastbeadclass.readcoords()[0,atom]) + ' ' + str(lastbeadclass.readcoords()[1,atom]) + ' ' + str(lastbeadclass.readcoords()[2,atom]) + '\n'
            output.write(entry.replace('e','E'))
            atom = atom + 1

    bead = bead + 1

print('')
print('The new ' + str(nbeads) + '-bead linear polymer has been written to: ' + str(outputfile)) 
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
