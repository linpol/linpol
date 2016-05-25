#!/bin/bash
# This file is part of LINPOL, Copyright 2013-2016 Michael F. Herbst
# (info@michael-herbst.com).
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

if [ "$1" ];then
	cd "$1"
fi

LOG=`ls -1t ????????-??????-log.txt | head -n1`

if [ -z "$LOG" ]; then
	echo "No log file found"
	exit 1
fi

overview() {
	echo "Overview of input params and chosen values:"
	TMP=`<$LOG  grep -A 6 "Reading atomic coordinates from xyz file"`
	echo "$TMP" | grep -B1 "beta is" | head -n1 | sed "s/^[[:space:]]*/   /g" 	#this is water XXmer
	BETA=`echo "$TMP" | grep "beta is" | head -n1 | awk 'BEGIN {FS = "beta is"}; {print $2}' | sed "s/.?//g"`
	TMP=`< $LOG  grep -A 3 "PARALLELIZATION: ALLOCATION OF BEADS ACROSS PROCESSES"`
	N=`echo "$TMP" | grep "There are" | awk '{print $3}'`
	printf "    %7s%15s       %10s%15s\n" "BETA = " "$BETA" "N       = " "$N"  | sed "s/^[[:space:]]*/   /g"

	if < $LOG  grep -q "^CALCULATING OPTIMUM GEOMETRY"; then
		TMP=`< $LOG  grep -A 6 "Running the L-BFGS minimization code to find instanton"`
		EPS=`echo "$TMP" | grep "EPS for optimization set to" | awk '{print $6}'`
		GRADEPS=`echo "$TMP" | grep "GRADEPS for bowman code set to" | awk '{print $7}'`
		printf "    %7s%15s       %10s%15s\n" "EPS  = " "$EPS"  "GRADEPS = " "$GRADEPS"  | sed "s/^[[:space:]]*/   /g"
	fi
	<$LOG  grep "^Using PES:" | head -n1 | sed  "s/^[[:space:]]*Using PES:/   PES =       /"
}

iteration() {
	if < $LOG grep -q "^Minimization finished"; then
		#Minimisation finished
		SUCC="without success"
		< $LOG  grep -q "*** L-BFGS CONVERGED SUCCESSFULLY ***" && SUCC="successfully"
		echo
		echo "Iteration converged $SUCC"
		echo "   Last iteration: `echo "$TMP" | tail -n1 | awk '{print $2}'`"
	else
		#not finished
		TMP=`cat $LOG | grep -E "^Iter"`
		echo
		echo "Iteration progression:"
		echo "   Last iterations of lowest gradient:"
		echo "$TMP" | grep  g | tail -n5 | sed "s/^/     /g"
		echo
		echo "   Last iterations of lowest function value:"
		echo "$TMP" | grep  f | tail -n5 | sed "s/^/     /g"
		echo
		echo "Iteration not finished"
		echo "   Current iteration:"
		echo "$TMP" | tail -n1 | sed "s/^/     /g"
	fi
}

propertycalc() {
	if < "$LOG" grep -q "Skink:"; then
		echo
		echo "Propery calculation results:"
		< "$LOG" grep "Bowman potential:" | sed "s/^[[:space:]]*/   /g"
		< "$LOG" grep "Skink:" | sed "s/^[[:space:]]*/   /g"
	fi
}

hessiancalc() {
	if < $LOG grep -q "^Diagonalizing Hessian matrix"; then
		if < $LOG grep -q "evalues ="; then
			echo
			echo "Hessian done"
		else
			echo
			echo "Currently diagonalising Hessian Matrix."
		fi
	else
		echo
		echo "Currently calculating Hessian Matrix:"
		< "$LOG" grep "Entry of Hessian:" | tail -n1 | sed "s/^[[:space:]]*/   /g"
	fi
}

splitting() {
	SPLITTING=`< "$LOG" grep -A8 "SPLITTING CALCUALTION" | sed "s/^[[:space:]]*/   /g"`
	echo "$SPLITTING" | grep  "log(det) ="| sed "s/^[[:space:]]*/   /g"
	echo "$SPLITTING" | grep  "phi ="| sed "s/^[[:space:]]*/   /g"
	echo "$SPLITTING" | grep  "h ="| sed "s/^[[:space:]]*/   /g"
}

###############################################################################
###############################################################################

overview
if < $LOG  grep -q "^CALCULATING OPTIMUM GEOMETRY"; then
	#a minimisation was done in this calculation
	iteration
fi
propertycalc

if < "$LOG"  grep -q "^Calculating Hessian matrix"; then
	#Hessian calculation was done in this calc
	hessiancalc
fi

if < "$LOG"  grep -q "SPLITTING CALCUALTION"; then
	splitting
fi

if < "$LOG" grep -q "Code complete"; then
	echo 
	<$LOG  grep -A6 "Runtime (per process) was:" 
fi
exit 0
