#!/usr/bin/python
# Filename: linpolmodule.py
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

from numpy import zeros

class datafile:
    def __init__(self,filetoopen):
        self.filetoopen = filetoopen
        self.fileobject = open(str(filetoopen) , 'r')
    def readlines( self ):
        self.fileobject.seek(0)
        linecount = 0
        for line in self.fileobject:
            linecount = linecount + 1
        self.fileobject.seek(0)
        self.nlines = int(linecount)
        return self.nlines
    def writelines( self ):
        self.readlines()
        print('Number of lines in ' + str(self.filetoopen) + ': ' + str(self.nlines))
    def readatoms( self ): 
        self.fileobject.seek(0)
        firstline = self.fileobject.readline()
        try:
            self.natoms = int(firstline)
        except ValueError:
            print('--------------------------------------------------------------------------')
            print('Error: Number of atoms is not an integer!')
            print('Check first line of ' + str(self.filetoopen))
            print('--------------------------------------------------------------------------')
            sys.exit('***SCRIPT ENDING PREMATURELY: Input file failure***')
        self.fileobject.seek(0)
        return self.natoms
    def writeatoms( self ):
        self.readatoms()
        print('Number of atoms in ' + str(self.filetoopen) + ': ' + str(self.natoms))
    def readnbeads( self ):
        self.readatoms()
        self.readlines()
        self.nbeads = self.nlines/(self.natoms+2)
        return self.nbeads
    def readbeta( self ):
        self.fileobject.seek(0)
        self.fileobject.readline()
        secondline = self.fileobject.readline()
        try:
            self.beta = float(secondline)
        except ValueError:
            print('--------------------------------------------------------------------------')
            print('Warning: Beta is not an floating point number or integer!')
            print('Check second line of ' + str(self.filetoopen))
            print('Script will continue, but please proceed with caution')
            print('--------------------------------------------------------------------------')
            #sys.exit('***SCRIPT ENDING PREMATURELY: Input file failure***')
            self.beta = secondline[:-1]		# To remove new line following
        self.fileobject.seek(0)
        return self.beta
    def writebeta( self ):
        self.readbeta()
        print('Beta in ' + str(self.filetoopen) + ': ' + str(self.beta))
    def readatomtags( self ):
        self.readatoms()
        self.atomtags = []
        self.fileobject.seek(0)
        self.fileobject.readline() # Skipping no. of atoms
        self.fileobject.readline()  # Skipping beta
        atom = 0
        while atom < self.natoms:
            string = self.fileobject.readline()
            string2 = string.split()
            self.atomtags.append(string2[0])
            atom = atom + 1
        return self.atomtags
    def readcoords( self ):
        self.readatoms()
        self.x = zeros((3,self.natoms))
        self.fileobject.seek(0)
        self.fileobject.readline() # Skipping no. of atoms
        self.fileobject.readline()  # Skipping beta
        atom = 0
        while atom < self.natoms:
            string = self.fileobject.readline()
            string2 = string.split()
            self.x[0,atom] = string2[1]
            self.x[1,atom] = string2[2]
            self.x[2,atom] = string2[3]
            atom = atom + 1
        return self.x
    def readall( self ):
        self.readatoms()
        self.x = zeros((4,self.natoms))
        self.fileobject.seek(0)
        self.fileobject.readline() # Skipping no. of atoms
        self.fileobject.readline()  # Skipping beta
        atom = 0
        while atom < self.natoms:
            string = self.fileobject.readline()
            string2 = string.split()
            if string2[0][0] == 'H':
                self.x[0,atom] = 0
            elif string2[0][0] == 'O':
                self.x[0,atom] = 1
            else:
                print('Unknown atom found (not hydrogen or oxygen)')
                sys.exit('***SCRIPT ENDING PREMATURELY: Input file failure***')
            self.x[1,atom] = string2[1]
            self.x[2,atom] = string2[2]
            self.x[3,atom] = string2[3]
            atom = atom + 1
    def readlinpol( self ):
        self.readatoms()
        self.readlines()
        self.readnbeads()
        self.x = zeros((4,self.natoms,self.nbeads))
        self.fileobject.seek(0)
        bead = 0
        while bead < self.nbeads:
            self.fileobject.readline() # Skipping no. of atoms
            self.fileobject.readline()  # Skipping beta
            atom = 0
            while atom < self.natoms:
                string = self.fileobject.readline()
                string2 = string.split()
                if string2[0][0] == 'H':
                    self.x[0,atom,bead] = 0
                elif string2[0][0] == 'O':
                    self.x[0,atom,bead] = 1
                else:
                    print('Unknown atom found (not hydrogen or oxygen)')
                    sys.exit('***SCRIPT ENDING PREMATURELY: Input file failure***')
                self.x[1,atom,bead] = string2[1]
                self.x[2,atom,bead] = string2[2]
                self.x[3,atom,bead] = string2[3]
                atom = atom + 1
            bead = bead + 1
        return self.x
