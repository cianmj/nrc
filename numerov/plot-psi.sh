#!/bin/sh
base="$1"
l="$2"
n="$3"
fils=$(echo $base-func-$l-$n.out|tr ' ' '\12'|sed -e 's/^/"/' -e 's/$/" using 1:($2\/$1), /' | tr '\12' ' '|sed -e 's/, *$//')
( tee xxx.gpl | gnuplot -persist ) <<EOI
ps=0
if (ps) set terminal postscript colour solid "Helvetica" 18
if (ps) set output "$base-func-$l-$n.ps"
set grid
plot $fils
EOI

