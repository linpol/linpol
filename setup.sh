#!/bin/bash
# This file is part of LINPOL, Copyright 2016 Michael F. Herbst
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

# Source this script in order to setup the environment for a linpol calculation.
_thisdir=$(dirname ${BASH_SOURCE[0]})

# Path to this directory (absolute):
export LINPOL="$(readlink -f "${_thisdir}")"

# Path to scripts (and linpolexe script, which sets up shop and calls the
# compiled linpolexe in the LINPOL dir)
export PATH="$PATH:$LINPOL/scripts"

# The default PES:
if [ -z "$LINPOL_DEFAULT_PES" ]; then
	export LINPOL_DEFAULT_PES="TTM3-F"
fi

# Disable scratch:
# The script linpolexe wraps around the actual linpolexe executable. 
# It automatically creates a unique scratch directory and performs 
# the calculation there. If this is not desired set:
# export LINPOL_DISABLE_SCRATCH="yes"	

unset _thisdir
