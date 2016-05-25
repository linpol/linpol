#!/bin/bash
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

#Settings & Defaults:
INPUT_GUESS="input_built.xyz"	#default file to be searched for in this directory and all subdirectories as bead input geometry
MINBETA=10000			#default minimum beta
MAXBETA=80000			#default maximum beta
MAXN=2048			#Maximum number of beads
MODE="double_beads"		#double_beads: starting guess for calculation (beta, N) is doubled bead geometry of (beta, N/2)
				#copy_geometry: starting guess for calculation (beta, N) is geometry of (beta/2,N)
EPS="1E-4"			#convergence tightness for L-BFGS
NPROC=1				#number of processes
[ -f /usr/bin/lscpu ] && NPROC=`lscpu | grep "^CPU(s):" | awk '{print $2}'`
TO="$USER"
EXTRA_FILTER='| grep -v "Iter:" | grep -v "Entry of Hessian:"'		#To suppress unwanted output from linpolexe in the exec.log (Note beginning | )
PES="$LINPOL_DEFAULT_PES"
STD_LINPOLEXE="linpolexe"
STD_CLUSTERDATA=""
STD_CHILDTO="$TO"		#where should mail from child mpigo processes sent to (only senceful option is "null", really)

#--------------------------------------------

prHelp() {
	cat << EOF
	`basename $0` [ Options ]

OPTIONS:
	-hall		Displays a longer help version.
	-h		Displays this short help.
	-t <TO>		Sends the email to the address <TO> instead of the
       			current user
	-nomail		Prevents child processes from sending mail.
	-p <No Procs>	Uses this number of processes (Default: Number of
       			Processors on current machine, or 1 if script can't 
			determine that. For your system the default
		       	would be $NPROC)
	-pes <PES>	Use this potential energy surface, default: $PES
	-data <DATA>	Read the cluster data from <DATA>
	-m <MODE>	either "double_beads" or "copy_geometry". 
			Defaults to the former. See section below for 
			more details.
	-minB <MINBETA>	minimum beta (default: $MINBETA)
	-maxB <MAXBETA> maximum beta (default: $MAXBETA)
	-maxN <MAXN>	maximum number of beads (default: $MAXN)
	-eps <EPS>	Standard convergence thightness for L-BFGS
			(default: $EPS)
	-i <INPUT>	Pattern for the input file (default: $INPUT_GUESS)
	-x <LINPOLEXE>  Use this linpolexe executable 
			(default: $STD_LINPOLEXE)

For more detailed help see "`basename $0` -hall"
EOF
}

prHelpAll() {
	cat << EOF
	`basename $0` [ Options ]

Sets up a LINPOL project to calculate the minimum energy instanton geometry for 
a series of BETA values and a series of linear polymer instantantons with 
variable size N.

The script determines the starting instanton guess by seeking the <INPUT> file
in the current directory and all subdirectories. By doubling the number of beads
and independantly the value of beta -- starting from <MINBETA> or the number of
beads in the <INPUT>, respectively -- a series of calculation is undertaken
until both <MAXBETA> and <MAXN> have been reached. 

For this purpose the script sets up a directory structure like
	<MINBETA>/<N>, <MINBETA>/<2*N>, ... <MINBETA>/<MAXN>
	<2*MINBETA>/<N>, <2*MINBETA>/<2*N>, ... <2*MINBETA>/<MAXN>
	...
	<MAXBETA>/<N>, <MAXBETA>/<2*N>, ... <MAXBETA>/<MAXN>
with one linpol calculation to be done in each directory. Parallel to this
process a bash shell script is written, which contains the sequence of commands
to be executed throughout the calculation. For more information about this script
see the section below the options.

OPTIONS:
	-hall		Displays this help.
	-h		Displays a short help.
	-t <TO>		Sends the email to the address <TO> instead of the
       			current user
	-nomail		Prevents child processes from sending mail.
	-p <No Procs>	Uses this number of processes (Default: Number of
       			processors on current machine, or 1 if script can't
		       	determine that. For your system the default
		       	would be $NPROC)
	-pes <PES>	Use this potential energy surface, default: $PES
	-data <DATA>	Read the cluster data from <DATA>
	-m <MODE>	either "double_beads" or "copy_geometry". 
			Defaults to the former. See section below for
		       	more details.
	-minB <MINBETA>	minimum beta (default: $MINBETA)
	-maxB <MAXBETA> maximum beta (default: $MAXBETA)
	-maxN <MAXN>	maximum number of beads (default: $MAXN)
	-eps <EPS>	Standard convergence thightness for L-BFGS
			(default: $EPS)
	-i <INPUT>	Pattern for the input file (default: $INPUT_GUESS)
			Current dir and all subdirs are parsed to find the file.
	-x <LINPOLEXE>  Use this linpolexe executable 
			(default: $STD_LINPOLEXE)

THE PROJECT MODES:
There are two principle modes a linpol project might run in. These are
controlled using the -m switch. "double_beads" means that the initial guess
geometry of a certain linear polymer (beta,N) with beta being the beta value for
the instanton and N being its number of beads is genereated by invoking the 
"0double.sh" script onto the minimum energy configuration of (beta,N/2), if 
available. "copy_geometry", however implies that (beta,N) gets its input geometry
from (beta/2,N). Hence in comparison we have the scenarios:

"double_beads"		<INPUT> == (beta,N)  --double-> (beta,2*N) --double->
			   |
			0change_beta
			  \/
			(2*beta,N)  -----------double-> (2*beta,2*N) --dble->


"copy_geometry"		<INPUT> == (beta,N)  --double-> (beta,2*N) --double->
			   |				    |
			0change_beta			0change_beta
			  \/				   \/
			(2*beta,N)			(2*beta,2*N)


THE $SKRIPT SCRIPT
The script executes the calculations in the individual subdirectories one after
another and checks if the results are sensible (using
0check_Linpol_Calculation.sh). If yes than the starting geometry for the next
calculation is produced and the next one executed. 
By default the Skript first executes the small N cases for all beta and then
moves on to the much longer calculations for higher beta. If the calculation
has stopped or did produce an error, it can be continued later. All folders
in which a "complete.txt" file is found are considered finished and will not
be touched.

NOTE -- The LINPOL_SCRIPTS variable:
Per default the $SKRIPT script fixes the LINPOL_SCRIPTS to the value it had
when `basename $0` was run.
The reason for this behavious is to make sure the correct underlying scripts
are called by the wrappers. Change the variable -before- you execute
$SKRIPT for the first time.
(In most cases it does, however, not make much of a difference what 
LINPOL_SCRIPTS is sassigned to)
EOF
}

calcMax() {
	#$1: start
	#$2: max
	#doubles $1 until max is reached or $1 is larger than max, returns this value
	RES=$1

	if [ "$2" -le "$1" ]; then
		echo $RES
		return 0
	fi

	while true;do
		RES=$[$RES*2]
		[ $RES -ge $2 ] && break
	done
	echo $RES
	return 0
}

check_configuration() {
	if [ -f "$SKRIPT" ];then
		echo "Script $SKRIPT exists."
	       	read -p "Do you want to continue and overwrite it? (y/n)  " RES
		[ -z "$RES" ] && RES="n"
		[ "$RES" != "y" ] && exit 1
	fi

	if type mail > /dev/null 2>&1;  then
		:
		#all is well
	else
		#ie no mail program
		echo "Error: There is no mail program installed on this system."
		echo "       Please ask your administrator to install one."
		exit 1
	fi

	if [ -z "$PES" ]; then
		echo "You did not select a proper potential energy surface."
		echo "Please provide one via -pes <PES>"
		exit 1
	fi

	if [ -z "$STD_CLUSTERDATA" ]; then
		echo "You did not provide proper cluster data."
		echo "Please provide a file via -data <DATA>"
		exit 1
	fi

	if [ -f "$INPUT_GUESS" ];then
		INPUT="$INPUT_GUESS"
	else
		INPUT=`find . -name "$INPUT_GUESS"`
		LC=`echo "$INPUT" | wc -l`
		if [ $LC -gt 1 ];then
			echo "Error: More than one file matched your input pattern \"$INPUT_GUESS\"."
			echo "       Please specify the total path of one of these file:"
			echo "$INPUT"
			exit 1
		fi
		if [ -z "$INPUT" ];then
			echo "Error: No file matched your input pattern \"$INPUT_GUESS\"."
			echo "       Please choose some other pattern or give the file "
			echo "       explicitly."
			exit 1
		fi
	fi

	#Determine number of beads:
	MINN=`cat "$INPUT" | awk '
		BEGIN { line1="";line2="";stage=0; c=0 }
		NR == 1 {line1=$0; next}
		NR == 2 {line2=$0;c++;next}; 
		$0 == line1 {stage=1;next};
		$0 == line2 && stage == 1 {c++;stage=0;next}	#increment if previous line was line1
		stage == 1 { stage=0 }				#previous line was line1, but this is not line2 => reset
		END{print c}'`

	if [ "$?" != "0" ]; then
		echo "Error opening file $INPUT"
		exit 1
	fi

	#for linpolexe the number of beads needs to be an integer multiple of the number of processors.
	WARN=0
	while true; do
		if [ "$[$MINN/$NPROC*$NPROC]" != "$MINN" ];then
			NPROC=$[NPROC-1]	#Note bash can only do integer arithmatic, ie $[3/2]=1
			WARN="1"
		else
			break
		fi
	done
	if [ "$WARN" = "1" ]; then
		echo "Warning: Number or processes was reduced to $NPROC to fit linpol requirements"
		echo
	fi

	if [ $MINN -gt $MAXN ]; then
		echo "Error: MAXN(= $MAXN) should be greater than the number of beads in the starting geometry(= $MINN)."
		exit 1
	fi
	MAXN=`calcMax $MINN $MAXN`

	if [ $MINBETA -gt $MAXBETA ]; then
		echo "Error: MAXBETA(= $MAXBETA) should be greater than MINBETA(= $MINBETA)."
		exit 1
	fi
	MAXBETA=`calcMax $MINBETA $MAXBETA`

	if [[ "$MODE" != "copy_geometry" && "$MODE" != "double_beads" ]];then
		echo "Error: $MODE is not a valid mode."
		exit 1
	fi
}

display_summary() {
	if [ "$STD_CHILDTO" != "null" ];then
		local MM="All mail gets sent to $TO."
	else
		local MM="A final summary gets sent to $TO."
	fi
	cat << EOF
The initial input geometry for (beta,N) == ($MINBETA,$MINN) is read from 
	$INPUT
The PES used is $PES and cluster data is read from $STD_CLUSTERDATA
The number of beads runs from $MINN to $MAXN.
Beta runs from $MINBETA to $MAXBETA.
The LBFGS convergence tightness is EPS = $EPS.
For the calculation we use $NPROC processes.
$MM

EOF
}

#--------------------------------------------
#build script

script_header() {
	MAXBETADFIGS=`echo -n $MAXBETA | wc -c`
	MAXNFIGS=`echo -n $MAXN | wc -c`

	( cat > $SKRIPT ) << EOF
#!/bin/bash

#-------------------------------------------------------
#Options:
export LINPOL_SCRIPTS="$LINPOL_SCRIPTS"
LINPOLEXE="$STD_LINPOLEXE"
PES="$PES"
CLUSTERDATA="$STD_CLUSTERDATA"
LOG="$PWD/exec_\`date +%Y%m%d-%H%M%S\`.log"
CHILDTO="$STD_CHILDTO"

#--------------------------------------------------------

read_config() {
	#read the config file in the current directory
	#and parse the args into the commandline args
	#passed to mpigo and linpolexe
	#the argline is saved in LINPOL_ARGS and INPUT

	local CONFIG=\`cat linpol_calc.config | sed "/^#/d"\`
	INPUT=\`echo "\$CONFIG" | awk 'BEGIN {FS="="}; \$1 ~ /^INPUT.?/ {print \$2 }' | tail -n1 | sed "s/[[:space:]]*//g"\`
	EPS=\`echo "\$CONFIG" | awk 'BEGIN {FS="="}; \$1 ~ /^EPS.?/ {print \$2 }' | tail -n1 | sed "s/[[:space:]]*//g"\`
	GRADEPS=\`echo "\$CONFIG" | awk 'BEGIN {FS="="}; \$1 ~ /^GRADEPS.?/ {print \$2 }' | tail -n1 | sed "s/[[:space:]]*//g"\`

	LINPOL_ARGS=""
	[ "\$GRADEPS" ] && LINPOL_ARGS="\$LINPOL_ARGS --gradeps \$GRADEPS"
	[ "\$EPS" ] && LINPOL_ARGS="\$LINPOL_ARGS --eps \$EPS"
	[ "\$INPUT" ] && LINPOL_ARGS="\$LINPOL_ARGS \$INPUT"
}

write_config() {
	#write the config file to the current directory
	#get the values from the corresponding vars

	( cat > linpol_calc.config ) << EOC
#Config file for the linpol calculation in this directory
#parsed by exec.sh
#lines starting with "#" are ignored
INPUT = \$INPUT
EOC
	[ "\$EPS" ] && echo "EPS = \$EPS" >> linpol_calc.config
	[ "\$GRADEPS" ] && echo "GRADEPS = \$GRADEPS" >> linpol_calc.config
}

send_log() {
	#\$1: some subject
	#Send log via mail.
	echo >> \$LOG
	echo >> \$LOG
	echo "Execution finished at \`date\`" >> \$LOG
	cat \$LOG | mail -s "[exec.sh]  \$1" $TO
}

#-------------------------------------------------------------

if [ "\$1" == "-h" ];then
	cat << EOD
The exec.sh script allows one of the following options
	--status	displays which folders have a finished calculation
	--stop		finishes the currently running calculation, then stops
	--clean		removes the complete.txt files from folders, unless
			they are marked to be kept. If a folder does not have
			a "complete.txt" it will be recalculated.
	--cleanall	removes all files from previous calculations from all
			folders but the ones which are marked as "keep".
	--tighten <LST>	Request tightening of convergence for the List given.
			Effectively the old optimised geometry is used as 
			input and the value for EPS is reduced by a factor
			of 10. (Only works with completed calculations 
			marked as "keep")
	--tightenall	request tightening for all calculations.
	--doonly <LIST>	only calculates the beta, No beads combinations
       			in the List <LIST>. Note: Dependancies are 
			not checked!
	--keep <LIST>	mark these calculations as ones to keep even for 
			invokation of "--clean" or "--cleanall". Effectively
			places a file "keep" in the folder. Also other calcula-
			tions won't overwrite the input.xyz file in these
			folders.
	--keepall	mark all calculations as "keep"
	--table		Runs straight "0extract_Linpol_table.sh --pdf" on calculation
			directory.
	--reinit <FILE> run "--cleanall" and also use <FILE> as the next initial
			guess for the first calculation ($MINBETA,$MINN).
	--vmdopt <N>	run vmdopt on each of the directores "BETA/N" for each
       			BETA. Good to check if results of a partcular number of
			beads have converged properly for various N.

LIST Format:
The Lists in the above commands have the format
	(B1,N1):(B2,N2):(B3,N3): ... :(BB,NN)
EOD
exit 0
fi

DOONLY=""

if [ -d "$PWD" ]; then
	cd "$PWD"
else
	echo "WARNING: could not set top directory properly:"
	echo "$PWD was not found"
	echo
fi

if [ "\$1" == "--status" ]; then
	COMPL=""
	NCOMPL=""
	for (( B=$MINBETA; B <= $MAXBETA  ;B=B*2 )); do
		for (( N=$MINN; N <= $MAXN; N=N*2 )); do
			cd "\$B/\$N/"
			read_config
			if [ -f "complete.txt" ];then
				extra=""
				LF=\`ls -1t ????????-??????-log.txt | head -n1\`
				if [ -f "keep" ];then
					extra="\$extra     marked as \\"keep\\""
				fi
				if grep "L-BFGS CONVERGED SUCCESSFULLY" "\$LF" > /dev/zero; then
					:
				else
					extra="\$extra     not converged properly"
				fi
				COMPL="\$COMPL:(\$B,\$N) to \$EPS \$extra"
			else
				NCOMPL="\$NCOMPL:(\$B,\$N) to \$EPS"
			fi
			cd ../..
		done
	done

	echo -n "The following calculations are completed:"
	echo "\$COMPL" | sed "s/:/\\n      /g"
	echo
	echo -n "The following are left to be done:"
	echo "\$NCOMPL" | sed "s/:/\\n      /g"
	exit 0
elif [ "\$1" == "--clean" ]; then
	for (( B=$MINBETA; B <= $MAXBETA  ;B=B*2 )); do
		for (( N=$MINN; N <= $MAXN; N=N*2 )); do
			if [ ! -f "\$B/\$N/keep" ]; then
				rm -f \$B/\$N/complete.txt
			fi
		done
	done
elif [ "\$1" == "--cleanall" ]; then
	for (( B=$MINBETA; B <= $MAXBETA  ;B=B*2 )); do
		for (( N=$MINN; N <= $MAXN; N=N*2 )); do
			if [ ! -f "\$B/\$N/keep" ]; then
				rm -f  \$B/\$N/*.txt \$B/\$N/*-*.xyz \$B/\$N/fort.*
			fi
		done
	done
	exit 0
elif [ "\$1" == "--tighten" ];then
	TOTIGHTEN=":\$2:"
	for (( B=$MINBETA; B <= $MAXBETA  ;B=B*2 )); do
		for (( N=$MINN; N <= $MAXN; N=N*2 )); do
			if echo "\$TOTIGHTEN" | grep ":(\$B,\$N):" > /dev/zero; then
				cd "\$B/\$N"
				if [[ -f "complete.txt" && -f "keep" ]]; then
					read_config #and fill variables

					PREFIX=\`ls -1 *-log.txt | head -n1\`
					PREFIX=\`basename \$PREFIX "-log.txt"\`
					mkdir "\$PREFIX"
					mv * "\$PREFIX" 2> /dev/zero

					cp "\$PREFIX/\$PREFIX-optimized.xyz" input.xyz
					echo "\$PREFIX/\$PREFIX-optimized.xyz" > input_from

					EPS=\`echo "\$EPS" | awk '{printf "%G", \$1/10}'\`

					write_config
					touch keep
				fi
				cd ../..
			fi
		done
	done
	exit 0
elif [ "\$1" == "--tightenall" ];then
	for (( B=$MINBETA; B <= $MAXBETA  ;B=B*2 )); do
		for (( N=$MINN; N <= $MAXN; N=N*2 )); do
			cd "\$B/\$N"
			if [[ -f "complete.txt" && -f "keep" ]]; then
				read_config #and fill variables

				PREFIX=\`ls -1 *-log.txt | head -n1\`
				PREFIX=\`basename \$PREFIX "-log.txt"\`
				mkdir "\$PREFIX"
				mv * "\$PREFIX" 2> /dev/zero

				cp "\$PREFIX/\$PREFIX-optimized.xyz" input.xyz
				echo "\$PREFIX/\$PREFIX-optimized.xyz" > input_from

				EPS=\`echo "\$EPS" | awk '{printf "%G", \$1/10}'\`
				echo "Tightening (\$B,\$N) to \$EPS requested"

				write_config
				echo keep
			fi
			cd ../..
		done
	done
	exit 0
elif [ "\$1" == "--keep" ]; then
	TOKEEP=":\$2:"
	for (( B=$MINBETA; B <= $MAXBETA  ;B=B*2 )); do
		for (( N=$MINN; N <= $MAXN; N=N*2 )); do
			if echo "\$TOKEEP" | grep ":(\$B,\$N):" > /dev/zero; then
				touch \$B/\$N/keep
			fi
		done
	done
	exit 0
elif [ "\$1" == "--keepall" ]; then
	for (( B=$MINBETA; B <= $MAXBETA  ;B=B*2 )); do
		for (( N=$MINN; N <= $MAXN; N=N*2 )); do
			touch \$B/\$N/keep
		done
	done
	exit 0
elif [ "\$1" == "--stop" ]; then
	touch stop_calc
	exit 0
elif [ "\$1" == "--table" ]; then
	0extract_Linpol_table.sh --pdf
	exit 0
elif [ "\$1" == "--doonly" ]; then
	DOONLY=":\$2:"
elif [ "\$1" == "--reinit" ]; then
	INPUT="\$2"
	if [ ! -f "\$INPUT" ]; then
		echo "Input file does not exist."
		exit 0
	fi

	./$SKRIPT --cleanall
	cp "\$INPUT" "$MINBETA/$MINN/input.xyz"
	0change_beta.sh "$MINBETA/$MINN/input.xyz" "$MINBETA"
	echo "\$INPUT (initial guess)" > "$MINBETA/$MINN/input_from"
	exit 0
elif [ "\$1" == "--vmdopt" ];then
	N="\$2"
	if [ ! -d "$MINBETA/\$N" ];then
		echo "Please provide a proper value for N"
      		exit 1
 	fi

	BETA=$MINBETA
	while [ \$BETA -le $MAXBETA ]; do
		cd "\$BETA/\$N"
		FILE=\`ls -1t *optimi?ed.xyz 2> /dev/zero | head -n1\`
		[ "\$FILE" ] && vmd \$FILE
		cd ../..
		BETA=\$[BETA*2]
	done
	exit 0		
elif [ "\$1" ]; then
	echo "Not a valid option: \$1."
	exit 1
fi

echo "Execution started at \`date\`" >> \$LOG
echo "ARGS: \$1 \$2" >> \$LOG
echo "LINPOLEXE: \$LINPOLEXE" >> \$LOG

#############################################
#############################################
#############################################

EOF
}

CONFIG_FILE_STANDARD_ACTION=""
script_drop_config() {
	local BETA=$1
	local N=$2
	#drops a config in the skript directory 

	if [ -f "$BETA/$N/linpol_calc.config" ]; then
		if [ "$CONFIG_FILE_STANDARD_ACTION" = "keepall" ]; then
			return
		fi
		if [ "$CONFIG_FILE_STANDARD_ACTION" != "overwriteall" ]; then
			echo
			echo "The file \"$BETA/$N/linpol_calc.config\" exists already."
			echo "Do you want to overwrite it with the standard config file?"
			echo "Type single character: y = Yes, n = No, a = All, x = none"
			read -n 1 -p "Choice   >  " RES
			echo

			case "$RES" in
				"x")	CONFIG_FILE_STANDARD_ACTION="keepall"
					return
					;;
				"a")	CONFIG_FILE_STANDARD_ACTION="overwriteall"
					;;
				"y")	CONFIG_FILE_STANDARD_ACTION=""
					;;
				*)	CONFIG_FILE_STANDARD_ACTION=""
					return
					;;
			esac
		fi
	fi	

	( cat > $BETA/$N/linpol_calc.config ) << EOF
#Config file for the linpol calculation in this directory
#parsed by exec.sh
#lines starting with "#" are ignored
INPUT = input.xyz
EPS = $EPS
#GRADEPS = 1E-4
EOF
}

script_element_begin() {
	local BETA=$1
	local N=$2
	#write the script part that does the calculation and checks the result.

	( cat >> $SKRIPT ) << EOF
if [ -f "stop_calc" ];then
	rm stop_calc
	echo "User requested STOP."
	echo "A stop was invoked." >> \$LOG
	send_log "User STOP request!"
	exit 0
fi

cd $BETA/$N

DOIT=0
if [ "\$DOONLY" ];then
	echo "\$DOONLY" | grep ":($BETA,$N):" > /dev/zero
	DOIT=\$?
fi
if [[ "\$DOIT" == "0" &&  ! -f "complete.txt" ]];then	#will be ended later!
	echo >> \$LOG
	echo "#######################################################################" >> \$LOG
	echo "#######################################################################" >> \$LOG
	echo "#######################################################################" >> \$LOG
	
	read_config

	#Not done yet & doonly requires us to do it
	if [ ! -f "\$INPUT" ]; then
		echo "Could not find a file \"\$INPUT\" in $BETA/$N" >> \$LOG
		echo "Aborting" >> \$LOG
		send_log "$BETA/$N: Couldn't find \$INPUT"
		exit 1
	fi

	echo "Invoke mpigo for $BETA/$N" >> \$LOG
	echo -n "Calculating ($BETA,$N)  "
	mpigo.sh -t \$CHILDTO -p $NPROC -x "\$LINPOLEXE" -- --data \$CLUSTERDATA --pes \$PES \$LINPOL_ARGS $EXTRA_FILTER &>> \$LOG
	if [ "\$?" != "0" ]; then
		echo "...  ERROR"
		echo "Error running mpigo.sh -x \"\$LINPOLEXE\"" >> \$LOG
		echo "Aborting" >> \$LOG
		send_log "$BETA/$N: mpigo/linpolexe failed"
		exit 1
	fi
	echo "...  done"

	echo "Checking results for $BETA/$N" >> \$LOG
	0check_Linpol_Calculation.sh "$PWD/$BETA/$N" >> \$LOG
	if [ "\$?" != "0" ]; then
		echo "Calculation not ok!" >> \$LOG
		echo "Aborting" >> \$LOG
		send_log "$BETA/$N: results check failed"
		exit 1
	fi
EOF
#if still open!!
}

script_element_end() {
	local BETA=$1
	local N=$2
	( cat >> $SKRIPT ) << EOF
fi					#if ends here !
cd ../..


EOF
}

script_element_double() {
	#write the script part that doubles the beads of the current geometry and saves the result in the appropriate dir.
	local BETA=$1
	local N=$2
	local twoN=$[$N*2]

	( cat >> $SKRIPT ) << EOF

	if [ -f "../$twoN/keep" ];then
		echo "Do not double beads from $N to $twoN, since keep file exists in $BETA/$twoN"  >> \$LOG
		echo "Skip doubling beads $N to $twoN, since keepfile"
	else
		OPT=\`ls -1t ????????-??????-optimized.xyz | head -n1\`
		0double.sh "\$OPT" "../$twoN/input.xyz"
		if [ "\$?" != "0" ]; then
			echo "Error doubling beads from $N to $twoN" >> \$LOG
			echo "Aborting" >> \$LOG
			send_log "$BETA/$N: error doubling beads"
			exit 1
		fi
		echo "$BETA/$N (double opt. geometry)" > "../$twoN/input_from"
	fi
EOF
}

script_element_copygeo() {
	#write the script part that copies the current bead geometry to the higher beta dir.
	local BETA=$1
	local N=$2
	local twoBETA=$[$BETA*2]

	( cat >> $SKRIPT ) << EOF

	if [ -f "../../$twoBETA/$N/keep" ];then
		echo "Do not copy geometry from $BETA/$N to $twoBETA/$N, since keep file exists in $twoBETA/$N" >> \$LOG
		echo "Skip copying geometry $BETA/$N to $twoBETA/$N, since keepfile"
	else
		OPT=\`ls -1t ????????-??????-optimized.xyz | head -n1\`
		cp "\$OPT" "../../$twoBETA/$N/input.xyz"
		if [ "\$?" != "0" ]; then
			echo "Error copying bead geometry from $BETA/$N to $twoBETA/$N" >> \$LOG
			echo "Aborting" >> \$LOG
			send_log "$BETA/$N: error copying bead geometry"
			exit 1
		fi
		0change_beta.sh "../../$twoBETA/$N/input.xyz" "$twoBETA"
		echo "$BETA/$N (opt. geometry)" > "../../$twoBETA/$N/input_from"
	fi
EOF
}

script_footer() {
	( cat >> $SKRIPT ) << EOF

echo >> \$LOG
echo "#######################################################################" >> \$LOG
echo "#######################################################################" >> \$LOG
echo "#######################################################################" >> \$LOG

send_log "ALL DONE"
EOF
}

########################################################################
########################################################################
########################################################################
SKRIPT="exec.sh"
INPUT=
MINN=

while [ "$1" ]; do
	case "$1" in
		"-hall")	prHelpAll
				exit 0
				;;
		"-h")	prHelp
			exit 0
			;;
		"-t")	shift
			TO="$1"
			;;
		"-nomail")
			STD_CHILDTO="null"
			;;
		"-p")	shift
			NPROC="$1"
			;;
		"-m")	shift
			MODE="$1"
			;;
		"-data") shift
			STD_CLUSTERDATA="$1"
			;;
		"-pes") shift
			PES="$1"
			;;
		"-minB") shift
			MINBETA="$1"
			;;
		"-maxB") shift
			MAXBETA="$1"
			;;
		"-maxN") shift
			MAXN="$1"
			;;
		"-eps")	shift
			EPS="$1"
			;;
		"-i") 	shift
			INPUT_GUESS="$1"
			;;
		"-x")   shift
			STD_LINPOLEXE="$1"
			;;
		*)	echo "Unknown argument $1"
			exit 1
			;;
	esac
	shift
done
#-----------------------------------------

check_configuration
display_summary
read -p "Press enter to start the creation of directories and scripts."

#-----------------------------------------

#Start script header:
script_header

N=$MINN
while [ $N -le $MAXN ];do
	BETA=$MINBETA
	while [ $BETA -le $MAXBETA ]; do
		mkdir -p $BETA/$N
		if [ "$?" != "0" ]; then
			echo "Error creating dir $BETA/$N"
			exit 1
		fi
		#----------------------------------
		script_element_begin "$BETA" "$N"
		script_drop_config "$BETA" "$N"

		if [ "$MODE" == "double_beads" ]; then
			#propagation along N
			if [ $[N*2] -le $MAXN ]; then
				script_element_double "$BETA" "$N"
			fi
			if [ $N == $MINN ];then
				if [ $[BETA*2] -le $MAXBETA  ]; then
					script_element_copygeo "$BETA" "$N"
				fi
			fi
		elif [ "$MODE" == "copy_geometry" ]; then
			#propagation along beta
			if [ $[BETA*2] -le $MAXBETA  ];then
				script_element_copygeo "$BETA" "$N"
			fi
			if [ $BETA == $MINBETA ]; then
				if [ $[N*2] -le $MAXN ]; then
					script_element_double "$BETA" "$N"
				fi
			fi
		fi
		script_element_end "$BETA" "$N"
		#----------------------------------
		BETA=$[BETA*2]
	done
	N=$[N*2]
done

if [ -f "$MINBETA/$MINN/input.xyz" ]; then
	echo
	echo "These already exists an initial guess input file in $MINBETA/$MINN/"
	read -p "Are you sure you want to overwrite it with the file $INPUT ? (y/n)  " RES
	[ -z "$RES" ] && RES="n"
else
	RES="y"
fi

if [ "$RES" == "y" ]; then
	cp "$INPUT" "$MINBETA/$MINN/input.xyz"
	0change_beta.sh "$MINBETA/$MINN/input.xyz" "$MINBETA"
	echo "$INPUT (initial guess)" > "$MINBETA/$MINN/input_from"
fi

script_footer
chmod 755 $SKRIPT

exit 0
