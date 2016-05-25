# This file is part of LINPOL, Copyright (C) Adam Reid, Jeremy Richardson 
# and Michael F. Herbst 2011-2016.
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
#
# Research conducted using LINPOL or any derivative thereof should cite 
# the following PhD thesis:
#
# A. A. Reid, "Quantum Tunnelling Splittings in Water Clusters, 
# from Ring-Polymer Instanton Theory", University of Cambridge (2014)
#
#########################################

#
# Compiler selection:
#
# Gfortran:
FC = mpif90 -O2 -fdefault-real-8
#
# Intel fortran:
#FC = mpiifort -O2 -ip -ipo -r8 -assume bscc -shared-intel

#
# Lapack selection:
#
# Plain reference lapack
LAPACK= -llapack
#
# AMD Core Math Library
# LAPACK = -lacml
#
# Threaded Intel MKL:
#LAPACK = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -liomp5
#
# Sequencial MKL:
#LAPACK = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

#
# Directories and files:
#
MBPOLDIR := src/mbpol
MBPOL_LIBS := $(MBPOLDIR)/mbpol/libmbpol.a
SRCDIR := src/NEB_01 src/linpol src/TTM3-F
SCRIPTSDIR=scripts
EXTRAFILES = README.md LICENCE example data
RELEASE=linpolexe_`date +%Y%m%d-%H%M%S`

#---------------

# Set libraries and fortran flags
LIBS := $(MBPOL_LIBS) ${LAPACK} -lstdc++
FFLAGS := 

#vpath tells make where to look for source files.
vpath %.f90 $(SRCDIR)
vpath %.f $(SRCDIR)

OBJ := 	math.o smear.o ttm3f_mod.o mnasa_mod.o mnasa.o phys_const.o\
	ttm3f.o init_ttm3f.o \
	pes_interface.o sort.o clusterdatamod.o parameters.o prep.o funcs.o \
	pots.o grads.o hess.o opt_events.o \
	type_module.o potential_optim_module.o lbfgs_module.o \
	opt_interface.o opt_selected.o args.o min.o linpolexe.o

# Targets:
.PHONY : clean release

%.o : %.f90
	$(FC) -c $(FFLAGS) $< 
%.o : %.f
	$(FC) -c $(FFLAGS) $< 

all : mbpol linpolexe

mbpol: 
	cd $(MBPOLDIR)/mbpol && $(MAKE)

linpolexe: $(OBJ) 
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS)

clean:
	rm -f *.o *.mod linpolexe *-log.txt; \
		cd $(MBPOLDIR)/mbpol && make clean

release:
	tar -cvzf $(RELEASE).tar.gz src $(SCRIPTSDIR) $(EXTRAFILES) \
		--exclude="*~"
