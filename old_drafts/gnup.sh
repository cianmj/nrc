#!/bin/sh
echo Creating PGN files...
cp density.txt dtemp.txt
cp densitybil.txt dbiltemp.txt
cp potential.txt ptemp.txt
max=$1
for ((i=0;$i<=$max;i=$i+1)) ; do
j=`echo $i | awk '{printf "%05d", $1}'`
echo $j
cat > temp.gpl <<-EOI
	set pointsize 0.2
	set xlabel "Position (au)"
	set autoscale
# Serguei
#       set xrange [-1.5:4.5]
#       set yrange [-0.0045:-0.002]
# Quartic Anharmonic
#        set xrange [-1.5:3.0]
#        set yrange [0.0:5.0]
# Harmonic
#       set xrange [-2.0:5.0]
#       set yrange [0.0:1.5]
# Morse
#       set xrange [0.5:2.5]
#       set yrange [0.0:0.05]
# Eckart
#       set xrange [-4.0:14.0]
#       set yrange [0.02:0.04]
# Hard Barrier
        set xrange [2.0:5.0]
        set yrange [0.0:1.1]
#
        set data style lines
	set output "images/$j.png"
	set terminal png
        plot "dtemp.txt" index $1, \
             "dbiltemp.txt" index $1, \
             "ptemp.txt"
EOI
gnuplot temp.gpl
done
rm ptemp.txt dtemp.txt dbiltemp.txt 
echo Finished.
