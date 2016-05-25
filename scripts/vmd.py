# This file is part of LINPOL, Copyright 2013-2016 Jeremy Richardson
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
"""Contains routines for inputing and outputing data in VMD readable format

Most commonly used functions are xyz and load.

"""

import numpy

#############################
# VMD Graphical representation settings:
# Method = CPK
# Material = Diffuse
# sphere resolution = 21
#
#############################

def getpdbline(atom, ID, x, charge):
	chargeval = abs(charge)
	if charge < 0: chargesym = '-'
	else: chargesym = '+'
	return "HETATM%5i %4s              %8.3f%8.3f%8.3f                      %2s%1i%1s" % (ID, atom, x[0], x[1], x[2], atom, chargeval, chargesym)

########################

# xyz format
#	  first line: number of atoms
#	  2nd line  : molecule name (can be blank)
#	  each atom with coords on new line
#	  multiple time steps

def getxyzline(atom, x):
	if x.size == 3:
		return "%s %.18g %.18g %.18g" % (atom, x[0], x[1], x[2])
	elif x.size == 2:
		return "%s %.18g %.18g 0.0" % (atom, x[0], x[1])
	elif x.size == 1:
		return "%s %.18g 0.0 0.0" % (atom, x)
	else:
		print "in vmd.getxyzline, x.shape =", x.shape
		exit(1)

def xyz(xt, atomlist, name, extra=None):
	"""Write coordinates in xyz format.
	xt - coordinates
		last dimension of xt must be 1, 2 or 3
		the middle dimensions fit to atomlist
		if present, the first dimension is time
	atomlist - nummpy.array of atomic symbols
	name - filename to be written
	extra - string to be placed on comment line (line 2 of file)
	"""
	xt = numpy.asarray(xt)
	atomlist = numpy.asarray(atomlist)

	lines = []

	if xt.ndim - atomlist.ndim == 1: # no time dimension
		lines.append(str(numpy.prod(xt.shape[0:-1])))
		if extra:
			lines.append(extra)
		else:
			lines.append(name)
		jshape = xt.shape[0:-1]
		for j in numpy.ndindex(jshape):
			j2 = tuple(numpy.where(numpy.array(atomlist.shape)==1, 0, j))
			lines.append(getxyzline(atomlist[j2], xt[j]))

	elif xt.ndim - atomlist.ndim == 2: # time dimension
		for t in range(xt.shape[0]):
			lines.append(str(numpy.prod(xt.shape[1:-1])))
			if extra:
				lines.append(extra)
			else:
				lines.append("t = %i" % t)
			jshape = xt.shape[1:-1]
			for j in numpy.ndindex(jshape):
				j2 = tuple(numpy.where(numpy.array(atomlist.shape)==1, 0, j))
				lines.append(getxyzline(atomlist[j2], xt[t][j]))

	else:
		print "ERROR in vmd.xyz: xt and atomlist have incompatible shapes"
		print "shape of xt is", xt.shape
		print "shape of atomlist is", atomlist.shape
		exit(1)

	with open(name, 'w') as f:
		for line in lines:
			f.write(line+"\n")

########################

def load(name, getatomlist=False):
	"""Read coordinates of xyz file."""
#	return numpy.loadtxt(name, skiprows=2, usecols=(1,2,3))
	with open(name, 'r') as f:
		data = f.readlines()
		N = int(data[0]) # number of atoms
		data = numpy.array(data)
		data.shape = (-1,N+2)	#hence each row is a bead; column1 = N, column2 = beta, column3 = atom1 symbol, 4 = atom1 xcoord ...
		data = data[:,2:]
		atomlist = numpy.empty(N, numpy.character)
		for j in range(N):
			atomlist[j] = data[0,j].split()[0]	#Assumes atom label is one character only!!
		x = numpy.empty((data.shape[0],N,3))
		for i in range(data.shape[0]):
			for j in range(N):
				x[i,j] = numpy.array(data[i,j].split()[1:], numpy.float)
	if x.shape[0] == 1:
		x = x[0]
	if getatomlist:
		return x, atomlist
	else:
		return x

########################

if __name__ == '__main__':
	
	if 1:
		x, atomlist = load("../../Instanton/H2OnHex/Bowman/dimer/trans/TS.xyz", True)
#		x = load("../../Instanton/H2OnHex/Bowman/dimer/trans/H2O+H2O/inst.32.350.xyz")
		print atomlist
		print x

	if 0: # xyz
		print getxyzline('H', numpy.array([1254.35, -0.00001, 1e5]))

	if 0: # pdb
		lines = []
		lines.append(getpdbline('H', 93, [0.1,0.2,0.3], 1))
		lines.append(getpdbline('H', 95, [0.1,0.2,1.3], 1))
		lines.append(getpdbline('O', 3353, [-0.1419, 2.2, 0.3], -2))
		print "12345678901234567890123456789012345678901234567890123456789012345678901234567890"
		for line in lines:
			print line
		with open('test.pdb', 'w') as f:
			for line in lines:
				f.write(line+"\n")
