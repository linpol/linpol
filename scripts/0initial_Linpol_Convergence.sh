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

START=
END=
HALF=		#optional halfway geometry
OPTSTART=0	#should starting geometry be optimised
OPTEND=0	#should end geometry be optimised
EPS="1E-5"	#epsilon for the instanton; for the geometries it is always EPS/100
EPSSTEPWISE=0	#reach EPS gradually when converging instanton
NPROC=`cat /proc/cpuinfo | grep processor | wc -l`		#number of processes for mpigo
DATA=""		#data file to use; if empty write data from one of the convergences of start, end
PES="$LINPOL_DEFAULT_PES"	#PES to use
TO=$USER	#To for various mpigo tasks
BETA=10000	#Beta to use for instanton
N=16		#Standard N to use
MAXITER=10000
FORCE="0"	#force overwrite

#------------------------------------------------------------------------

STARTFOLDER="optimised_startpoint"
ENDFOLDER="optimised_endpoint"
INSTANTFOLDER="initial_convergence"
FOLDERDATAFILE="beaddata.dat"
EPSSTEPNUMBER=4  #number of calculations to do in a stepwise calculation
		 #starts at <EPS>*10^<NUMBER> and proceeds to tighten in powers of 10

#------------------------------------------------------------------------

prHelp() {
	cat << EOF
	`basename $0` [ Options ] <XYZ Start> <XYZ End>

Creates a Folder "$INSTANTFOLDER" in which an instanton guess is 
set up and converged. The guess contains
	<N>/2 start and 
	<N>/2 end beads
if NO HALFWAY bead guess is given and 
	<N>/2-<N>/16 start,
	<N>/8        halfway and
	<N>/2-<N>/16 end beads
if a HALFWAY bead guess is given.

OPTIONS:
	-h		Displays this short help.

	--optstart	optimise the start bead before using it.
			Will create a folder "$STARTFOLDER".
	--optend	optimise the end bead before using it.
			Will create a folder "$ENDFOLDER".
	--half <FILE>	Also include a guess for the halfway geometry
			in the guess for the instanton.

	-b <BETA>	Beta for the instanton guess. (DEFAULT: $BETA)
	-n <N>		Number of beads for the instanton guess (DEFAULT: $N)
	--eps <EPS>	Convergence tightness for the instanton guess.
			(default: $EPS)
			For the bead geometries always <EPS>/100 is taken.
	--stepwise	Tighten the instanton guess stepwise
			Will produce a sequence of $EPSSTEPNUMBER calculations where
			slowly but steadily EPS is tightened. 

	--maxiter <IT>  Use a maximum number of <IT> iterations for each
       			calculation.

	--pes <PES>	Use this potential energy surface, default: $PES
	--data <DATA>	Read the cluster data from <DATA>. If not given and 
			either start or end is optimised we use the data 
			generated from that optimisation. Start is preferred
			if both are meant to be optimised. If none are and 
			--data is missing and error is printed.

	-t <TO>		Sends the email to the address <TO> instead of the
       			current user
	--nomail	Prevents child processes from sending mail.
	-p <No Procs>	Uses this number of processes (Default: Number of
       			Processors on current machine. For your system 
			the default would be $NPROC)
	--force		Force overwriting some data
EOF
}

###########################################

check_configuration() {
	if [[ -z "$START" || -z "$END" ]]; then
		echo "Please provide both a valid starting and a valid end bead geometry"
		exit 1
	fi
	
	#number of atoms in start bead:
	local STARTA=`< $START grep -xEc "^[[:space:]]*[[a-zA-Z]]?[[:space:]]*-*[0-9]?\.[0-9]?.*[[:space:]]*-*[0-9]?\.[0-9]?.*[[:space:]]*-*[0-9]?\.[0-9]?.*$"` 
	#number of atoms in end bead:
	local ENDA=`< $END grep -xEc "^[[:space:]]*[[a-zA-Z]]?[[:space:]]*-*[0-9]?\.[0-9]?.*[[:space:]]*-*[0-9]?\.[0-9]?.*[[:space:]]*-*[0-9]?\.[0-9]?.*$"` 

	if [ "$ENDA" != "$STARTA" ]; then
		echo "Start and End bead need to have the same number of atoms"
		exit 1
	fi

	if [ "$HALF" ]; then
		local HALFA=`< $HALF grep -xEc "^[[:space:]]*[[a-zA-Z]]?[[:space:]]*-*[0-9]?\.[0-9]?.*[[:space:]]*-*[0-9]?\.[0-9]?.*[[:space:]]*-*[0-9]?\.[0-9]?.*$"` 
		if [ "$HALFA" != "$ENDA" ]; then
			echo "Halfway bead needs to have the same number of atoms as start and end"
			exit 1
		fi
	fi
	
	if [ -z "$PES" ]; then
		echo "Please provide an potential energy surface to use."
		echo "Either use --pes or set the environment variable \$LINPOL_DEFAULT_PES"
		exit 1
	fi

	if [ -z "$DATA" ]; then
		if [ "$OPTSTART" = "1" ];then
			STARTWRITEDATA=1
			DATA="$PWD/$STARTFOLDER/$FOLDERDATAFILE"
		elif [ "$OPTEND" = "1" ];then
			ENDWRITEDATA=1
			DATA="$PWD/$ENDFOLDER/$FOLDERDATAFILE"
		else
			echo "No data file specified and no optimisation will be run."
			echo "Please do either provide a file or optimise end/start"
			exit 1
		fi
	else
		if [ -f "$DATA" ]; then
			#make sure we have an absolute path the the data file
			DATA="$PWD/$DATA"
		fi
	fi

	if [ "$FORCE" != "1" ]; then
		if [[ "$OPTSTART" == "1" && -d "$STARTFOLDER" ]]; then
			echo "Optimisation required, but $STARTFOLDER exists."
			echo "Please move folder to prevent data to be overwritten or --force"
			exit 1
		fi
		if [[ "$OPTEND" == "1" && -d "$ENDFOLDER" ]]; then
			echo "Optimisation required, but $ENDFOLDER exists."
			echo "Please move folder to prevent data to be overwritten or --force"
			exit 1
		fi
		if [ -d "$INSTANTFOLDER" ]; then
			echo "Instanton guess folder exists."
			echo "Please move folder to prevent data to be overwritten or --force"
			exit 1
		fi
	fi
}

display_summary() {
	echo "#################################"
	echo "#-- Summary of input provided --#"
	echo "#################################"
	echo

	local TMP="gets optimised to eps=`echo "$EPS" | awk '{printf "%G", $1/100}'`"
	[ "$OPTSTART" != "1" ] && TMP="does NOT get optimised"
	[ "$STARTWRITEDATA" == "1" ] && TMP="$TMP and cluster data will be written to $DATA" 
	echo "Startfile     $START $TMP"

	local TMP="gets optimised to eps=`echo "$EPS" | awk '{printf "%G", $1/100}'`"
	[ "$OPTEND" != "1" ] && TMP="does NOT get optimised"
	[ "$ENDWRITEDATA" == "1" ] && TMP="$TMP and cluster data will be written to $DATA" 
	echo "Endfile       $END $TMP"
	[ "$HALF" ] && echo "Halfwayfile   $HALF"

	echo "The instanton guess will be optimised according to:"
	echo "  BETA    =   $BETA"
	echo "  N       =   $N"
	echo "  PES     =   $PES"
	local GRADUALLY=; [ "$EPSSTEPWISE" = "1" ] && GRADUALLY="(gradually)"
	echo "  EPS     =   $EPS  $GRADUALLY"
	echo "  DATA    =   $DATA"
	echo
	read -p "Press enter to start the optimisation"
	echo
	echo
}

###############################################

optbead() {
	#$1: file with guess
	#$2: folder to create and carry out optimisation
	#$3: shall we write data: 1 for yes
	
	echo "#################################"
	echo "#-- Optimising bead $1"
	echo "#################################"
	echo

	local OEPS=`echo "$EPS" | awk '{printf "%G", $1/100}'`
	local ARGS=" --maxiter $MAXITER --pes $PES --eps $OEPS guess.xyz"
	if [ "$3" = "1" ]; then
		ARGS="--writedata $DATA $ARGS"
	else
		#we just need the geometry
		ARGS="--nohess $ARGS"
	fi

	mkdir -p "$2"
	cp $1 "$2/guess.xyz"
	# go into subshell:
	(
		cd "$2"
		eval "linpolexe $ARGS"

		if [ "$?" != "0" ]; then
			echo "Linpolexe returned with an error"
		  	exit 1
		fi
	)
	[ "$?" != "0" ] && exit 1 
	echo
	echo
}

build_instantion() {
	NSTARTEND=$((N/2)) 
	NHALF=$((N/8))
	[ "$HALF" ] && NSTARTEND=$((NSTARTEND-N/16))

	mkdir -p "$INSTANTFOLDER"

	cp "$START" "$INSTANTFOLDER/start.xyz"
	cp "$END" "$INSTANTFOLDER/end.xyz"
	[ "$HALF" ] && cp "$HALF" "$INSTANTFOLDER/half.xyz"

	rm -f "$INSTANTFOLDER/instanton_guess.xyz"


	#add starting geometry at the beginning:
	for (( i=1; i <= $NSTARTEND; i++)); do
		cat "$INSTANTFOLDER/start.xyz" >> "$INSTANTFOLDER/instanton_guess.xyz"
	done

	if [ "$HALF" ];then
		for (( i=1; i <= NHALF; i++)); do
			cat "$INSTANTFOLDER/half.xyz" >> "$INSTANTFOLDER/instanton_guess.xyz"
		done
		NSTARTEND=$((NHALF+NSTARTEND))
	fi

	#add end geometry at the end:
	for ((i=$NSTARTEND+1;i <= $N;i++)); do
		cat "$INSTANTFOLDER/end.xyz" >> "$INSTANTFOLDER/instanton_guess.xyz"
	done

	0change_beta.sh "$INSTANTFOLDER/instanton_guess.xyz" "$BETA"
}

opt_instanton() {
	#$1: Folder
	#$2: EPS to use
	#$3: Extra Args

	echo "#################################"
	echo "#-- Optimising instanton $1"
	echo "#################################"
	echo

	#go to subshell:
	(
		cd "$1"
		mpigo.sh -t $TO -p $NPROC -x linpolexe -- --maxiter $MAXITER --pes $PES --eps $2 --data "$DATA" $3 input.xyz
		if [ "$?" != "0" ]; then
			echo "Linpolexe returned with an error"
		  	exit 1
		fi

		0check_Linpol_Calculation.sh .
		if [ "$?" != "0" ]; then
			echo "Linpolexe Calculation result check failed"
		  	exit 1
		fi
	)
	[ "$?" != "0" ] && exit 1 
	echo
	echo
}

########################################################################
########################################################################
########################################################################

STARTWRITEDATA=0
ENDWRITEDATA=0

while [ "$1" ]; do
	case "$1" in
		"--optstart")
			OPTSTART=1
			;;
		"--optend")
			OPTEND=1
			;;
		"--stepwise")
			EPSSTEPWISE=1
			;;
		"--force")
			FORCE=1
			;;
		"-h")	prHelp
			exit 0
			;;
		"-t")	shift
			TO="$1"
			;;
		"--nomail")
			TO="null"
			;;
		"-p")	shift
			NPROC="$1"
			;;
		"-b")	shift
			BETA="$1"
			;;
		"-n")	shift
			N="$1"
			;;
		"--half")
			shift
			if [ -f "$1" ]; then
				HALF="$1"
			else
				echo "Invalid xyz file for halfway bead: $1"
				exit 1
			fi
			;;
		"--data") shift
			DATA="$1"
			;;
		"--maxiter") shift
			MAXITER="$1"
			;;
		"--pes") shift
			PES="$1"
			;;
		"--eps") shift
			EPS="$1"
			;;
		*)	if [ -f "$1" ]; then
				if [ "$START" ]; then
					END="$1"
				else
					START="$1"
				fi
			else
				echo "Unknown argument or invalid file $1"
				exit 1
			fi
			;;
	esac
	shift
done
#-----------------------------------------

check_configuration
display_summary

if [ "$OPTEND" -eq "1" ];then
	optbead "$END" "$ENDFOLDER" $ENDWRITEDATA
	END=`echo $ENDFOLDER/????????-??????-optimi?ed.xyz`
fi

if [ "$OPTSTART" -eq "1" ];then
	optbead "$START" "$STARTFOLDER" $STARTWRITEDATA
	START=`echo $STARTFOLDER/????????-??????-optimi?ed.xyz`
fi

build_instantion

cd "$INSTANTFOLDER"
#convert INSTANTFOLDER to absolute path:
INSTANTFOLDER="$PWD"

#---

[ "$EPSSTEPWISE" == "0" ] && EPSSTEPNUMBER=1

CEPS=`echo "$EPS" | awk -v "n=$EPSSTEPNUMBER" '{printf "%8.2E", $1*10^(n-1)}'`
PREV="instanton_guess.xyz"

#runs EPSSTEPNUMBER-1 times:
for (( i=1; i < $EPSSTEPNUMBER; i++)); do
	mkdir -p "$CEPS"
	cp "$PREV" "$CEPS/input.xyz"
	echo "../$PREV" > "$CEPS/input_from"
	opt_instanton "$CEPS" "$CEPS" "--nohess"

	PREV=`echo "$CEPS/"????????-??????-optimi?ed.xyz`
	CEPS=`echo "$CEPS" | awk '{printf "%8.2E", $1/10}'`
done

mkdir -p "$EPS"
cp "$PREV" "$EPS/input.xyz"
echo "../$PREV" > "$EPS/input_from"
opt_instanton "$EPS" "$EPS"
PREV=`echo "$EPS/"????????-??????-optimi?ed.xyz`


echo "Initial convergence procedure finished successfully."
echo "Find the converged instanton at"
echo "$INSTANTFOLDER/$PREV"
exit 0
