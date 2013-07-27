#!/bin/sh
if [ $1 = pos ]
 then
  num=2
elif [ $1 = mom ]
 then
  num=3
elif [ $1 = nrg ]
 then
  num=4
else
 echo Incorrect input
 exit 0
fi
cat > "$1".gpl <<-EOI
	set pointsize 0.2
	set xlabel "time (fs)"
	set autoscale
       set xrange [0.0:2500.0]
#       set yrange [0.0:1.5]
	set output "iplots/$1-$2.png"
	set terminal png
	plot "optimal.10" using 1:$num with lines, \
             "classical.11" using 1:$num with lines
#	pause -1
EOI
gnuplot "$1".gpl
cd iplots
xv $1-$2.png
cd ..
rm *.gpl
