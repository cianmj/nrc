#!/bin/sh
file="give a name, bozo"
[ -r "$1-eigenvalues.out" ] && file="$1-eigenvalues.out"
[ -r "$1.levels"          ] && file="$1.levels"
awk < $file \
'
BEGIN { 
   first = 1 ; 
   print "set grid"
   print "set samples 400"
#  print "set xrange [0:10]"
#  print "set xrange [10:20]"
#  print "set xrange [20:30]"
#  print "set xrange [30:44]"
   print "set xrange [0:45]"
   }
!/^#/ { 
   l = $1 ; n = $2 ; z = $3 ; e = $4 ; 
   if (first) {
      printf "plot " ;
      first = 0 ;
      }
   else {
      printf ", \\\n     "
      }
   printf "(x<%d||x>%d)?1/0:%5.3f notitle with lines %d", l, l+1, e, l+1 ;
   }
' | tee $file.gpl | gnuplot -persist
