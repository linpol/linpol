#!/bin/sh
BASE=$(basename "$0" ".pdf.sh")
MD="$BASE.md"
pandoc --smart --number-sections --standalone --toc -V "linkcolor:red" \
	-M geometry="scale=.8" -o "$BASE.pdf" "$@" "$BASE.md"
