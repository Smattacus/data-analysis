#!/usr/bin/sh
#Remakes the latex file given as an argument with bibtex.

bibname=`echo $1 | sed -e "s/.tex//g"`
dviname=`echo $1 | sed -e "s/.tex/.dvi/g"`
echo $bibname
echo $dviname
latex $1
bibtex $bibname
latex $1
latex $1
dvipdf $dviname
echo "dvipdf $dviname"
