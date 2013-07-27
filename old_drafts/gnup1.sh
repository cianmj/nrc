#!/bin/sh
cat > test.gpl <<-EOI
	set pointsize 0.2
	set xlabel "Position (au)"
	set autoscale
# Serguei
#       set xrange [-1.5:4.5]
#       set yrange [-0.0045:-0.002]
# Quartic Anharmonic
#	set xrange [-1.5:3.0]
#	set yrange [0.0:5.0]
# Harmonic
#	set xrange [-2.0:5.0]
#	set yrange [0.0:1.5]
# Morse
#	set xrange [0.5:2.5]
#	set yrange [0.0:0.05]
# Eckart
#	set xrange [-4.0:14.0]
#	set yrange [0.02:0.04]
# Hard Barrier
	set xrange [1.0:4.5]
	set yrange [0.0:1.1]
#
	set data style lines
	plot "density.txt" index $1, \
         "densitybil.txt" index $1, \
         "potential.txt"
	pause -1
EOI
gnuplot test.gpl
exit 0
