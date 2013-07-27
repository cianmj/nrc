#!/bin/sh
file="$1"
lval="$2"
gnuplot -persist <<EOI
set grid
plot "$file-score-$lval.out" using 1:(log(\$3)) with linespoints
EOI
