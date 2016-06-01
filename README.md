# LINPOL: LINear POLymer instanton code

This document serves as a small introduction to the ``linpol`` program.
If you find this work useful please consider citing the references given below.  

This document can also be viewed in pdf format by executing the script
[README.pdf.sh](README.pdf.sh).

## About LINPOL and its authors
LINPOL is a LINear POLymer instanton code for calculating tunnelling matrix
elements. It was developed at the University of Cambridge by Adam Reid and others
between 2011 and 2013.

LINPOL makes use of an amended form of the L-BFGS optimization algorithm
[Liu and J. Nocedal., Math. Program. 45, 503 (1989)].

Other contributors to the project have included:

- *Jeremy Richardson*: a source of general advice on the implementation of
instanton theory and a contributor of some pre-existing subroutines for
the calculation of the Hessian matrix, plots of tunnelling pathways, etc.

- *Michael F. Herbst*: contributor of structural enhancements to the code,
parallelization of the Hessian calculation routine and implementation of
routines to calculate phi and h.  Michael also created many of the useful
scripts which accompany the program.

- *Marko Cvitaš*: contributor of the amended L-BFGS algorithm used by LINPOL.

Other code (not written by Reid et al) included alongside LINPOL:

- ``src/mbpol``: MP-pol potential of V. Babin, G. R. Medders and F. Paesani

- ``src/TTM3-F``: TTM3-F potential of G. S. Fanourgakis and S. S. Xantheas

It should be noted that this code does not form part of LINPOL and it is
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. The copyright remains with the original authors as noted above or in the
corresponding source files.

## Compiling and setup
By default the ``Makefile`` is setup to compile LINPOL and the afforementioned
libraries with:

- The GNU compilers ``gfortran`` and ``g++``

- The default ``lapack`` library.

If other compilers (like the Intel compilers) or other LAPACK implementations
(Intel MKL, ACML) are desired, the Makefiles need to be amended manually.
See the Makefiles ``Makefile`` and ``src/mbpol/mbpol/Makefile`` for details.

Once the setup has been done properly, the code can be compiled with a simple
``make``.

### Setting up the environment
In order to use ``linpolexe`` from a terminal some environment variables, like
the location of linpolexe and a default potential energy surface, need to
be set up. For this purpose source the script ``setup.sh``. It chooses ``TTM3-F`` by
default. Once you have run
```
. /path/to/LINPOL/setup.sh
```
from a terminal you are good to go. More information on the environment
variables and how to use the ``linpolexe`` program specifically can be found
by running
```
linpolexe -hall
```
after the setup.

By default all calcultions of ``linpolexe`` are not done in the current working
directory but much rather in ``/tmp`` or ``/scratch/$USER``. If this behaviour
is not desired, it can be disabled by setting ``LINPOL_DISABLE_SCRATCH`` to
``yes``. See the file ``setup.sh`` for details.

## Getting started
LINPOL calculates $S_{kink}$ (the action) and $\phi$ (the ratio of determinants),
and thus $h$ (the tunnelling matrix element), for a particular tunnelling
pathway between two versions of a molecule (named version 1 and version 2 below).
To do  this, it requires the following:

- Environment variables to have been set (see above).

- The creation of a single-bead data file containing information about the
  two versions of the cluster which are being studied.

- Some pregenerated data is available in the ``data`` subfolder for ``TTM3-F``.

- To add single-bead data file for a cluster version the user must supply an
  ``input.xyz`` geometry file, which must then be optimized by LINPOL. The data
  is then written to file, which is specified by ``--writedata``:
  ```
  linpolexe --pes <PES> --writedata $LINPOL/data/<PES>/cluster.dat input.xyz
  ```
  The process is relatively quick and has not been parallelised.

Once these prerequisites are in place a calculation of a single instanton is
started as follows:
```
linpolexe --data <STRING> --pes <PES> --eps <EPS> --hesseps <HESSEPS> \
    <STARTING_GUESS_FOR_OPTIMISATION.xyx>
```

For more details about the command line flags which can be provided to
``linpolexe``, see the output of ``linpolexe -hall``.

## Converging a matrix element
What follows is an edited excerpt from A. A. Reid's PhD thesis (2014).

In order to calculate the tunnelling matrix element for the pathway between versions
1 and 2 of a molecule/cluster, the values of $S_{kink}$ and $phi$ must be numerically
converged to large $\beta$ (i.e. low temperature) and large $N$ (the number of beads).
In practice, a convenient and effective way of proceeding with a set of calculations
is as follows:

- Create a starting geometry of 64 beads in which the first 32 beads represent version 1,
and the last 32 beads represent version 2.

- Run the L-BFGS search on this geometry at $\beta \hbar$ = 40,000

- Use the resulting optimized geometry as a starting guess for the 64-bead calculations
at $\beta \hbar$ = 20,000 and $\beta \hbar$ = 80,000.

- Double the beads in the optimized 64-bead geometries at all three $\beta \hbar$ values
(use the useful ``0double.sh`` script) and use the resulting 128-bead geometries as starting
points for the next round of L-BFGS searches.

If the pathway in question is suitable for analysis with ring-polymer instanton theory,
the results of a process of this nature will resemble those in the tables below:

Table of values of $S_{kink}$:

$N / \beta \hbar$  |  20,000  |  40,000  |  80,000
------------------ |  ------- |  ------- |  -------
64                 |  23.02   |  22.36   |  12.30
128                |  23.15   |  23.02   |  23.80
256                |  23.18   |  23.15   |  23.02
512                |  23.19   |  23.18   |  23.15
1024               |  23.19   |  23.19   |  23.18
2048               |  23.19   |  23.19   |  23.19
4096               |  23.19   |  23.19   |  23.19

Table of values of $\Phi$:

$N / \beta \hbar$  |   20,000  |  40,000  | 80,000
------------------ |  -------- | -------- |  -------
64                 |   10.6    |  70.6    |  13671.4
128                |   15.0    |  10.6    |  2420.2
256                |   22.7    |  14.9    |  10.5
512                |   28.0    |  22.5    |  14.7
1024               |   30.1    |  27.8    |  22.3
2048               |   30.6    |  29.8    |  27.5
4096               |   30.8    |  30.4    |  29.5

A study of such tables allows us to check for convergence of both $S_{kink}$ and $\Phi$; once we
are happy that this has occured, the tunnelling matrix element is found from the converged
results ($\beta \hbar$ = 80,000 and $N$ = 4,096).

## Scripts to help with the creation of convergence tables.
In his master's thesis Michael F. Herbst wrote a bunch of shell scripts which aid
with setting up the calculations necessary for creating such tables.
For easy tab completion these files all start with a ``0``. Help and a much more
detailed description of their functionality is available by providing the ``-h``
flag. Try for example
```
0setup_Linpol_Project.sh -h
```

A short description of the most important scripts:

- ``0initial_Linpol_Convergence.sh``: Form an initial starting geometry for the
    tunneling calculation. Only guesses for the start and end geometry need to
    be available. The script takes care of optimising these in the desired
    potential as well as forming the inital geometry.

- ``0setup_Linpol_Project.sh``: Setup a project with all calculations required
    to form tables for $S_{kink}$ and $\Phi$. This script will generate an
    ``exec.sh`` script and a number of directories for the calculations. These
    can then be orchestrated using ``exec.sh``. Detailed help is available with
    ``0setup_Linpol_Project.sh -hall``

- ``0extract_Linpol_table.sh`` extract $S_{kink}$ and $\Phi$ data from the
    calculations, make a tex file with the table and build a pdf.

For more details about the computational strategies that led to the
development of these scripts, see his master's thesis, which is available for
download from
[http://blog.mfhs.eu/2013/08/10/master-thesis-tunnelling-in-water-clusters/]
(http://blog.mfhs.eu/2013/08/10/master-thesis-tunnelling-in-water-clusters/).

## Example input files
These are contained in the ``example`` subdirectory.
Try optimizing them using LINPOL. They can be viewed in VMD.

## Including further potentials:
Further potentials can be included by editing the source file
``src/linpol/pes_interface.f90`` and linking the potential into ``linpolexe``
during compilation.

## LINPOL references
Research conducted using LINPOL or any derivative thereof should cite the
following PhD and master's thesis:

- A. A. Reid, *Quantum Tunnelling Splittings in Water Clusters, from
    Ring-Polymer Instanton Theory*, University of Cambridge (2014)
- M. F. Herbst, *Instanton calculations of tunnelling in water clusters*,
    University of Cambridge (2013)

For general help regarding ring-polymer instanton theory, the central point
of contact is Professor Stuart Althorpe, ``sca10@cam.ac.uk``
