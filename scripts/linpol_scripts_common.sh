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

##################################################################################
###############
#-- OpenMPI --#
###############

mpi_is_root() {
	#is this the root process
	if [[ -z "$OMPI_COMM_WORLD_RANK" || "$OMPI_COMM_WORLD_RANK" == "0" ]]; then
		return 0
	else
		return 1
	fi
}

#################################################################################
###############
#-- Scratch --#
###############
#Manages the scratch for calculations

get_scratch_top() {
	if [ -d "/scratch/$USER" ];then
		echo "/scratch/$USER"
		return 0
	fi
	echo "/tmp/${USER}_scratch"
}

goto_scratch() {
	#All of this is run by EACH mpi process!!

	#setup
	export SCRATCH_OLDDIR=`pwd`

	if [[ $SCRATCH_OLDDIR =~ "/tmp" ]]; then
		#we are in TMP => need no scratch
		unset SCRATCH_OLDDIR
		return
	fi

	#Set LINPOL_DISABLE_SCRATCH="yes" in bashrc to disable switching to scratch
	if [ "$LINPOL_DISABLE_SCRATCH" == "yes" ]; then
		unset SCRATCH_OLDDIR
		return
	fi

	if mpi_is_root; then
		SCRATCH_DATE=`date +%Y%m%d-%H%M%S.%N`			#when was scratch created
		export SCRATCH_DIR="`get_scratch_top`/$SCRATCH_DATE"	#the scratchdir
		mkdir -p "$SCRATCH_DIR"

		#copy & create files
		echo "$SCRATCH_OLDDIR" > "$SCRATCH_DIR/.scratch_clone_from"
		cp -ru * "$SCRATCH_DIR"

		#This has to be the last thing (tells other processes that we are ready):
		echo "$SCRATCH_DIR" > "$SCRATCH_OLDDIR/.scratch_clone_at"
	else
		#make sure that root has copied all the files and everything is set up ready to go:
		local C=0
		while [ ! -f "$SCRATCH_OLDDIR/.scratch_clone_at" ]; do
			sleep 1
			if [ $[++C] -gt 100 ];then
				echo "FATAL ERROR: It took the root process more than 100 seconds toset up the scratch dir."
				echo "I assume an error occurred there and will exit now."
				exit 1
			fi
		done
		export SCRATCH_DIR=`cat "$SCRATCH_OLDDIR/.scratch_clone_at"`
	fi
	#-----------------------------------
	cd "$SCRATCH_DIR"
}

leave_scratch() {
	[ -z "$SCRATCH_OLDDIR" ] && return #scratch not set up!

	if mpi_is_root; then
		rm "$SCRATCH_DIR/.scratch_clone_from"
		cp -ru * "$SCRATCH_OLDDIR"
		rm -rf "$SCRATCH_DIR"
		rm "$SCRATCH_OLDDIR/.scratch_clone_at"
	fi
	cd "$SCRATCH_OLDDIR"
}

##################################################################################
##################################################################################
##################################################################################

#--------------------------------
#test prerequisites
if [ -z "$LINPOL" ]; then
	echo "Please set environment variable \$LINPOL to point to the top directory containing the linpol executable."
	exit 1
else
	export PYTHONPATH="$PYTHONPATH:$LINPOL/scripts"
fi

