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

PATTERN=
OUT="movie.avi"

prhelp() {
	cat << EOF
	`basename $0` [ -h ] [ -o <FILE> ] <FILE_PATTERN>

Make a movie using ffmpeg discarding
  - N/10  frames at the beginning
  - N/5   frames in the middle
  - N/10  frames at the end
where N is the number of files matching the pattern. 
The file pattern contains a construct "%0nd" 
where n is the number of digits, hence "%05d" would match 
00000, 00001, 00002, ..., 99999. The requirement is, however, 
that the files are numbered consecutively.

Options:
-o <FILE>	name of the final movie file. (Default: $OUT)
EOF
}

#-------------------------------------------------------

while [ "$1" ];  do
	case "$1" in
		"-h")	prhelp
			exit 1
			;;
		"-o")	shift
			OUT="$1"
			;;
		*)	if echo "$1" | grep -Eq "%0.d"; then
				PATTERN="$1"
			else
				echo "Unknown switch or invalid Pattern: $1"
				exit 1
			fi
			;;
	esac
	shift
done

#-------------------------------------------------------

if [ -z "$PATTERN" ]; then
	echo "Please provide a valid pattern"
	exit 1
fi

#-------------------------------------------------------

function number_to_file() {
	#$1: a number
	#This function adds trailing zeros to the number and substitutes it into pattern.

	local NBR="$ZEROS$1"
	NBR="${NBR: -$DIGITS}"

	echo "$PATTERN" | sed "s/%0${DIGITS}d/$NBR/"
}

#-------------------------------------------------------

#Determine number of digits from pattern and build zero string:
DIGITS=`echo "$PATTERN" | grep -o "%0.d"`; DIGITS=${DIGITS:2:1}
ZEROS=`echo "$DIGITS" | awk '{for (i=0;i<$1;i++) s=s "0"}; END{print s}'`


#determine first and last:
FIRST=0
LAST=
echo -n "finding the start and end  ...  "
while true; do
	[ -f `number_to_file $FIRST` ] && break
	((FIRST=FIRST+1))
done
LAST=$((FIRST+1))
while true; do
	[ ! -f `number_to_file $LAST` ] && break
	((LAST=LAST+1))
done
((LAST=LAST-1))
echo "done"

N=$((LAST - FIRST + 1))
FROM1=$((N/10))
TO1=$((N/2-N/5))
FROM2=$((TO1+N/5))
TO2=$((N-N/10))

#-----------------------------

echo -n "Moving files               ...  "
mkdir -p left_out

C=0
while (( C < FROM1 )); do
	mv `number_to_file $C` left_out/
	((C=C+1))
done

#offset by which we shift files:
OFFS=$((FROM2 - TO1 - 1))
C=0
while (( (MOVEFROM=FROM2+C) <= TO2 )); do
	MOVETO=$((MOVEFROM-OFFS))
	mkdir -p left_out
	FLE=`number_to_file $MOVETO`
	[ -e "$FLE" ] && mv `number_to_file $MOVETO` left_out/
	mv `number_to_file $MOVEFROM` `number_to_file $MOVETO`

	((C=C+1))
done

((C=TO2+1))
while (( C < TO2-1 )); do
	mv `number_to_file $C` left_out/
	((C=C+1))
done
echo "done"
echo

#-----------------------------

FROM1="$ZEROS$FROM1"
FROM1="${FROM1: -$DIGITS}"
ffmpeg -f image2 -start_number $FROM1 -i "$PATTERN" -b:v 500k "$OUT"

