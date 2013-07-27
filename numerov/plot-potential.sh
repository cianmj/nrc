#!/bin/sh
#
#  Plot Intra-cage potential
#
system="$1"
#
[ ! -r "$system.levels"    ] && { echo "$system.levels" not found ; exit 1 ; }
[ ! -r "$system.potential" ] && { echo "$system.potential" not found ; exit 1 ; }
#
#
#
awk </dev/null \
'
function pow(x,n,  v,i){
  v = 1 ; for (i=1;i<=n;i++){ v *= x ; }
  return v ;
  }
'"$(cat "$system.potential")"'
END {
  for (r=0;r<=rmax;r+=rs) {
     printf " %15.5f  %25.15f\n", r, pot(r) ;
     }
  }
' rmax="$(awk '/RMAX/{print $2}' $system.levels)" rs=0.01 > "$system.potential.dat"
#
#  Plot radial distribution
#
( tee "$system.potential.gpl" | gnuplot -persist ) <<EOI
set grid
set xlabel "r, Angstrom"
set ylabel "v, kcal/mol"
plot "$system.potential.dat" title "v(r)" with lines 1
EOI
