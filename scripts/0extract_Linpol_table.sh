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

#Settings & Default:
OUT="table.tex"
CHECK="y"
OVERWRITE="y"
PDF="n"

prHelp() {
	cat << EOF
	`basename $0` [ Options ] [ <FILE> ]

Searches the subdirectoris for linpol calculation directories with finished
calculations and perfors some checks in order to determine if calculation
has finished properly:
	- is "complete.txt" present
	- run 0check_Linpol_Calculation.sh on directory
	- is "keep" keepfile present (see 0setup_Linpol_Project.sh)
		(implies that results have been checked)

If the tests are passed the scrpt greps for the splitting parameters 
(phi, Skink and h) in the log file and parses these into three LaTeX tables
saved as Skink_<FILE>, Phi_<FILE> and h_<FILE>.

OPTIONS:
	-h		Displays this help.
	--nocheck	Do not perform any proper checks. If logfile exists
			try to parse it.
	--pdf		Also create a pdf file which displays the table
EOF
}

parse_dirs() {
	#parses the directories and saves the importand data in the temporary file $TMP.
	#The different N and BETA values found are put into the file $TMP.2 in the format:
	#          $NS#-#-#-#-#-#$BETAS
	echo "Parsing directories"
	find . -name "????????-??????-log.txt" | sort -g -k 2  --field-separator=/ --stable | sort -g -k 3 --field-separator=/ --stable | { 
		while read file; do
			#file contains names of logfiles
			DIR=`dirname "$file"`

			if [ "$CHECK" == "y" ]; then
				0check_Linpol_Calculation.sh "$DIR"
				if [ "$?" != "0" ];then
					echo "   ... skipped \"$file\" (check failed)"
					continue
				fi

				if [ ! -f "$DIR/keep" ]; then
					echo "   ... skipped \"$file\" (no \"keep\" file)"
					continue
				fi
			fi

			N=`cat "$file" | grep -A5  "Reading values from file" | awk ' $0 ~ /^These lines were split into/ {print $6}'`
			BETA=`cat "$file" | grep -A5  "Reading atomic coordinates from xyz file" | awk '$1 == "beta" && $2 == "is" {print $3}'`
			SKINK=`cat "$file" | grep -A5  "PROPERTY CALCUALTION" | grep  "Skink:" | cut -d ":" -f 2`
			PHI=`cat "$file" | grep -A8  "SPLITTING CALCUALTION" | grep  "phi =" | cut -d "=" -f 2`
			H=`cat "$file" | grep -A8  "SPLITTING CALCUALTION" | grep  "h =" | cut -d "=" -f 2`
			
			echo -e "$BETA\t$N\t$SKINK\t$PHI\t$H" >> $TMP

			if echo "$NS" | grep -v " $N " > /dev/zero ; then
				#not yet added:
				NS="$NS $N "
			fi
			if echo "$BETAS" | grep -v " $BETA " > /dev/null; then
				#not yet added:
				BETAS="$BETAS $BETA "
			fi

			echo "   ... added \"$file\""
		done
		echo "$NS#-#-#-#-#-#$BETAS" > "$TMP.2"
	} # end of subshell after find pipe
}
#--------------------------------------
table_header() {
	#$1:	NS
	#prints the table header

	HEAD=""
	COLS="r|"
	for N in $NS; do
		COLS="${COLS}c"
		HEAD="$HEAD & $N"
	done
	HEAD="$HEAD \\\\"


cat << EOF
\\begin{tabular}{$COLS}
	$HEAD
	\hline
EOF
}

table_row(){
	#$1:	ROW
	echo "	$1 \\\\" | sed "s/#-#-#-#-#-#/ \& /g"
}

table_footer() {
cat << EOF
\\end{tabular}
EOF
}

wrapper_PDF() {
	#$1: OUT
	#prints a small wrapper pdf to STDOUT
cat << EOF
\\documentclass[landscape,a4paper]{article}
\usepackage{amsmath}
\usepackage{lscape}
\usepackage[tight]{units}
\\begin{document}

\setlength {\\textheight}{160mm}
\setlength {\\textwidth}{250mm}
\setlength {\\oddsidemargin}{0mm}
\setlength {\\topmargin}{-15mm}

\\section*{\$S_{\\text{kink}}\$ Table (in units of \$ \hbar \$)}
\\begin{center}
\\input{Skink_$1}
\\end{center}

\\section*{\$\Phi\$ Table (in a.u. of time)}
\\begin{center}
\\input{Phi_$1}
\\end{center}


\\section*{\$h_{\lambda \mu}\$ Table (kink path weight; in Hartree)}
\\begin{center}
\\input{h_$1}
\\end{center}

Note:
\\[ 1 E_h \\equiv \unit[219474.63]{cm^{-1}} \\equiv \unit[627.5095]{kcal/mol} \\equiv \unit[6.579\,683\,920\,729(33) \\cdot 10^{15}]{Hz}  \\]

\\end{document}
EOF
}

############################################################################
############################################################################
############################################################################
TMP=`mktemp`

while [ "$1" ]; do
	case "$1" in
		"-h")	prHelp
			exit 0
			;;
		"--nocheck")
			CHECK="n"
			;;
		"--nooverwrite")
			OVERWRITE="n"
			;;
		"--pdf")
			PDF="y"
			;;
		*)	if [ -f "$1" ]; then
				OUT="$1"
			else
				echo "Unknown argument and not a valid File  $1"
				exit 1
			fi
			;;
	esac
	shift
done

parse_dirs
if [ "$?" != "0" ]; then
	exit 1
fi

#$TMP.2:	$NS#-#-#-#-#-#$BETAS
NS=`cat "$TMP.2" | awk 'BEGIN {FS="#-#-#-#-#-#"}; {print $1}'`
BETAS=`cat "$TMP.2" | awk 'BEGIN {FS="#-#-#-#-#-#"}; {print $2}'`


#sort the file:
sort "$TMP" > "$TMP.2"
mv "$TMP.2" "$TMP"

#Tables:
table_header "$NS" > Skink_$OUT
table_header "$NS" > Phi_$OUT
table_header "$NS" > h_$OUT
for BETA in $BETAS; do
	ROW="$BETA"	#separated by #-#-#-#-#-#
	ROWPHI="$BETA"	#separated by #-#-#-#-#-#
	ROWH="$BETA"	#separated by #-#-#-#-#-#
	for N in $NS; do
		# Jeremy: 
		# I recommend printing the Skink table with only two decimal places as the 
		# error in the 2nd decimal place is only exp(0.01)-1=1%. This is the sort 
		# of accuracy we're looking for and you won't get anything better because 
		# the numerical optimizer is never perfect.I recommend printing the Skink 
		# table with only two decimal places as the error in the 2nd decimal place 
		# is only exp(0.01)-1=1%. This is the sort of accuracy we're looking for 
		# and you won't get anything better because the numerical optimizer is 
		# never perfect.
		#=> use %6.2f (2 figs after decimal)
		SKINK=`cat $TMP | awk -v "n=$N" -v "b=$BETA" '$1 == b && $2 == n { printf("%6.2f",$3) }'`
		[ -z "$SKINK" ] && SKINK="--"

		# Jeremy:
		# The table of Phi should be printed to 2 significant figures as this
		# again gives about 1% errors.
		#=> %6.2G (2 sig. figs)
		PHI=`cat $TMP | awk -v "n=$N" -v "b=$BETA" '$1 == b && $2 == n { printf("%6.2G",$4) }'`
		[ -z "$PHI" ] && PHI="--"

		#For the output of H I used 2 signifficant figures.
		H=`cat $TMP | awk -v "n=$N" -v "b=$BETA" '$1 == b && $2 == n { printf("%6.2G",$5) }'`
		[ -z "$H" ] && H="--"

		ROW="${ROW}#-#-#-#-#-#${SKINK}"
		ROWPHI="${ROWPHI}#-#-#-#-#-#${PHI}"
		ROWH="${ROWH}#-#-#-#-#-#$H"
	done
	table_row "$ROW" >> Skink_$OUT
	table_row "$ROWPHI" >> Phi_$OUT
	table_row "$ROWH" >> h_$OUT
done
table_footer >> Skink_$OUT
table_footer >> Phi_$OUT
table_footer >> h_$OUT

rm $TMP
if [ "$PDF" == "y" ];then
	OLDDIR=`pwd`
	mkdir $TMP
	cp Skink_$OUT Phi_$OUT h_${OUT} $TMP
	cd $TMP
	PDFNAME="`basename "$OUT" ".tex"`.pdf"
	wrapper_PDF "$OUT" > $OUT
	pdflatex $OUT
	if [ "$?" == "0" ];then
		mv $PDFNAME $OLDDIR
	else
		echo 
		echo "pdflatex ended with an error!"
	fi
	cd $OLDDIR
	rm -r $TMP
fi

exit 0
