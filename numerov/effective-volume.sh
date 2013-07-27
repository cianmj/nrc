#!/bin/sh
#
#  Calculates classical vs. quantum effective volumes at a given 
#  temperature range.
#
system="$1"
t1="$2"
t2="$3"
ts="$4"
#
[ ! -r "$system.levels"    ] && { echo "$system.levels" not found ; exit 1 ; }
[ ! -r "$system.potential" ] && { echo "$system.potential" not found ; exit 1 ; }
#
#  Quantum part - calculate partition function Z, then multiply by the
#                 cube of thermal wave length.
#
awk < "$system.levels" \
'
BEGIN {
  gasR    = 0.001987 ; # kcal/mol/k
  lambdaA = 17.46 ; # Angstrom*AMU**0.5*K**0.5
  npt     = 2000 ;
  n       = 0 ;
  }
function pow(x,n,  v,i){
  v = 1 ; for (i=1;i<=n;i++){ v *= x ; }
  return v ;
  }
'"$(cat "$system.potential")"'
/^#\$MASS/ { m = $2 ; printf "# Using mass = %8.2f AMU\n", m ; }
/^#\$RMAX/ { rmax = $2 ; printf "# Using rmax = %8.4f Angstrom\n", rmax ; }
!/^#/ { 
  lval = $1 ; nroot = $2 ; energy = $4 ;
  n++ ;
  eps[n] = $4 ;
  wgt[n] = 2*$1 + 1 ;
  }
function quantumq(t,  lambda,i,z) {
  lambda = lambdaA / sqrt(m*t) ;
  z      = 0 ;
  for (i=1;i<=n;i++){
     z += wgt[i] * exp(-eps[i]/(gasR*t)) ;
     }
  return pow(lambda,3) * z ;
  }
function classicalq(t,   q,ipt,r,expv) { # Naive numerical integration
  q = 0 ;
  for (ipt=1;ipt<npt;ipt++) {  # Midpoints
    r = ipt*rmax/npt ;
    q   += r*r*exp(-pot(r)/(gasR*t)) ;
    }
  q += 0.5*rmax*rmax*exp(-pot(rmax)/(gasR*t)) ;
  q *= 4 * 3.1415926 ;
  return q*rmax/npt ;
  }
END {
  printf "# Volumes are in Angstrom**3\n" ;
  printf "# %8s  %15s  %15s\n", "T, K", "Classical", "Quantum"
  for (t=t1;t<=t2;t+=ts) {
     printf " %8.1f  %25.4f  %25.4f\n", t, classicalq(t), quantumq(t) ;
     }
  }
' t1="$t1" t2="$t2" ts="$ts" > "$system.veff.dat"
#
( tee "$system.veff.gpl" | gnuplot -persist ) <<EOI
set grid
set ylabel "q_{eff}, Angstrom^3"
set xlabel "T, Kelvin"
plot "$system.veff.dat" using 1:2 title "classical" with lines 1, \
     "$system.veff.dat" using 1:3 title "quantum"   with lines 2
EOI
#
( tee "$system.veffr.gpl" | gnuplot -persist ) <<EOI
set grid
set ylabel "(q_{classical}-q_{quantum})/q_{quantum}, percent"
set xlabel "T, Kelvin"
plot "$system.veff.dat" using 1:((\$2/\$3-1.)*100.) title "classical" with lines 1
EOI
