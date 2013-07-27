#!/bin/sh
base="$1"
l="$2"
n="$3"
fils=$(echo $base-func-$l-$n.out|tr ' ' '\12'|sed -e 's/^/"/' -e 's/$/" using 1:($2), /' | tr '\12' ' '|sed -e 's/, *$//')
gnuplot -persist <<EOI
set grid
plot $fils
EOI

