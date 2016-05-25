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
"""For converting ring polymers to normal mode representation

cyclic - for ring polymers
linear - for linear polymers

use forward and backward for speedy operations
transform is slower but gives normal modes in same form as given in the appendix of JCP 131 214106.

For a usage example, run this file as a script.

"""

import math, numpy
from scipy import fftpack

#############################################################

def backwardsbroadcast(a, x):
	"""append extra axes to a to be able to broadcast with x"""
	a.shape += tuple(numpy.ones(x.ndim-a.ndim, numpy.int))
	return a

#############################################################

class cyclic:
	def transform(self, x):
		shape = x.shape
		N = x.shape[0]
		x.shape = (N,-1)
		q = numpy.empty_like(x)
		q[0] = numpy.sum(x,0) / math.sqrt(N)
		for k in range(1,N/2):
			q[k] = math.sqrt(2./N) * numpy.sum(x * util.backwardsbroadcast(numpy.sin(2*math.pi*numpy.arange(N)*k/N), x), 0)
			q[-k] = math.sqrt(2./N) * numpy.sum(x * util.backwardsbroadcast(numpy.cos(2*math.pi*numpy.arange(N)*k/N), x), 0)
		if N % 2 == 0: # N is even
			q[N/2] = numpy.sum(x * util.backwardsbroadcast((-1)**numpy.arange(N), x), 0) / math.sqrt(N)
		x.shape = shape
		return q.reshape(shape)
		
	def forward(self, x):
		shape = x.shape
		N = x.shape[0]			#No of beads
		x.shape = (N,-1)
		q = numpy.empty_like(x)
		for j in range(x.shape[1]):	#normal mode analysis as fft
			q[:,j] = fftpack.rfft(x[:,j])
		q /= math.sqrt(N)
		if N % 2 == 0: # N is even
			q[1:-1] *= math.sqrt(2)		#all but the first and last one
		else:
			q[1:] *= math.sqrt(2)		#all but first one
		x.shape = shape			#bring x and q to the shape (No beads, No atoms in bead, 3) - 3 for xyz coords
		return q.reshape(shape)
		
	def backward(self, q):
		N = q.shape[0]
		x = q.copy()
		x.shape = (N,-1)
		x *= math.sqrt(N)
		if N % 2 == 0: # N is even
			x[1:-1] /= math.sqrt(2)
		else:
			x[1:] /= math.sqrt(2)
		for j in range(x.shape[1]):
			x[:,j] = fftpack.irfft(x[:,j])
		return x.reshape(q.shape)
	
	def matrix(self, N):
		X = numpy.identity(N)
		return self.forward(X)
	
#############################

class linear(cyclic):
	def transform(self, x):
		q = numpy.empty_like(x)
		N = x.shape[0]
		for k in range(N):
			q[k] = math.sqrt(2./(N+1)) * numpy.sum(x * util.backwardsbroadcast(numpy.sin(math.pi*numpy.arange(1,N+1)*(k+1)/(N+1)), x), 0)
		return q

	def forward(self, x):
		return self.transform(x)

	def backward(self, q):
		return self.transform(q)

#############################

if __name__ == "__main__":
	NM = cyclic()
	x = numpy.array([0, 0.55, 0.6, 0.4, -0.1, -0.6, -0.7, -0.33])
#	x = numpy.array([[-1,0.1], [0,0], [1,-0.1], [0,0]], numpy.float)
	print "x", x
	J = NM.matrix(len(x))
	print "transform"
	q = NM.transform(x)
	print q
#	print numpy.dot(J, x)
#	print q[0], math.sqrt(q[1]**2+q[-1]**2), math.sqrt(q[2]**2+q[-2]**2), math.sqrt(q[3]**2+q[-3]**2), q[4]
	print "forward"
	q = NM.forward(x)
	print q
#	q = numpy.concatenate((q, numpy.zeros_like(q))) * math.sqrt(2)
#	print q
#	print q[0], math.sqrt(q[1]**2+q[2]**2), math.sqrt(q[3]**2+q[4]**2), math.sqrt(q[5]**2+q[6]**2), q[7]
	print "backward"
	x = NM.backward(q)
	print x
#	print numpy.dot(J.T, q)
#	print numpy.mean(x)

	print "\nlinear"
	NM = linear()
#	x = numpy.array([0, 0.1, 0.24, 0.3, 0.5, 0.6, 0.7, 0.75])
	x = numpy.array([-1, 0, 1], numpy.float) + 1
	print "x", x
#	print "transform"
#	q = NM.transform(x)
#	print q
	print "forward"
	q = NM.forward(x)
	N = q.size
	q = numpy.concatenate((q, numpy.zeros(N))) * math.sqrt((2*N+1.0)/(N+1.0))
	print q
	print "backward"
	x = NM.backward(q)
	print x
	print numpy.mean(x)
