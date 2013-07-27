#!/bin/sh
#
#  Calculates classical vs. quantum radial distribution function
#  for a given temperature.
#
system="$1"
temp="$2"
#
[ ! -r "$system.levels"    ] && { echo "$system.levels" not found ; exit 1 ; }
[ ! -r "$system.potential" ] && { echo "$system.potential" not found ; exit 1 ; }
#
#  Quantum part - calculate Boltzman factors for all known levels,
#                 and construct corresponding ensemble average.
#
awk < "$system.levels" \
'
BEGIN {
  gasR      = 0.001987 ; # kcal/mol/k
  psiscl    = 1.374673 ; # normalization factor
  sumwgt    = 0 ;
  npoints   = 0 ;
  lmax      = -1 ;
  nmax      = -1 ;
  have_grid = 0 ;
  eps       = 3e-5 ;
  }
function abs(x) { return x<0?-x:x ; }
function loadFunction(l,n,  name,thisn,norm) {
  name = sprintf("%s/%s-func-%d-%d.out",base,base,l,n) ;
  thisn = 0 ; norm  = 0 ;
  while ( (getline < name) > 0 ){
     if (substr($0,1,1) != "#" ) {
        thisn = thisn + 1 ;
        if (!have_grid) {
           rval[thisn] = $1 ; pval[thisn] = 0 ; npoints = thisn ;
           }
        else {
           if (abs($1-rval[thisn])>eps) {
              printf "Grid mismatch in file %s\n", name > "/dev/stderr" ; exit(1) ;
              }
           }
        fv[thisn] = psiscl*$2 ; norm += fv[thisn] * fv[thisn] ;
        }
     }
  if (thisn==0) {
     printf "File %s does not exist or is empty\n", name ;
     exit(1) ;
     }
  if (thisn!=npoints) {
     printf "Number of grid points changed in %s\n", name ;
     exit(1) ;
     }
  have_grid = 1 ;
  norm *= (rval[npoints]-rval[1]) / npoints ;
  if (abs(norm-1.0)>0.01){
     printf "Function %s is normalized to %12.7f\n", name, norm > "/dev/stderr" ;
     }
  }
function addState(l,n,w,  i) {
  loadFunction(l,n)

  for (i=1;i<=npoints;i++) {
     pval[i] += w*fv[i]*fv[i]
     }
  }
!/^#/ { 
  lval = $1 ; nroot = $2 ; energy = $4 ;
  scale   = (2*lval+1)*exp(-energy/(gasR*temp)) ;
  if (lval >lmax) lmax = lval ;
  if (nroot>nmax) nmax = nroot ;
  wgt[lval,nroot] = scale ;
  sumwgt += scale ;
  addState(lval,nroot,scale) ;
  }
END {
  printf "## Weights of individual states (degeneracy-adjusted):\n" ;
  printf "##  L  N  Weight\n" ;
  for (lval=0;lval<=lmax;lval++) {
     for (n=0;n<=nmax;n++) {
        if (wgt[lval,n] > 0 ){
           printf "#  %3d %3d %20.13f\n", lval, n, wgt[lval,n]/sumwgt ;
           }
        }
     }
  rs = 0 ; ps = 0 ;
  printf "#%12s  %12s\n", "R, Angstrom", "Probability" ;
  for (i=1;i<=npoints;i++){
     printf " %12.6f  %12.6f\n", rval[i], pval[i]/sumwgt ;
     if (i<npoints) {
        ps += (pval[i]/sumwgt) * (rval[i+1]-rval[i]) ;
        rs += (rval[i+1]+rval[i]) * (pval[i]/sumwgt) * (rval[i+1]-rval[i]) / 2 ;
        }
     }
  printf "Quantum: total probability = %12.6f, <r> = %12.6f\n", ps, rs > "/dev/stderr" ;
  }
' temp="$temp" base="$system" > "$system.quantum.$temp.dat"
#
#  Classical part - calculate Boltzman weight factors for the same
#                   potential on the same grid
#
awk < "$system.quantum.$temp.dat" \
'
BEGIN {
  swgt = 0 ;
  gasR = 0.001987 ; # kcal/mol/k
  n    = 0 ; swgt = 0 ;
  }
function pow(x,n,  v,i){
  v = 1 ; for (i=1;i<=n;i++){ v *= x ; }
  return v ;
  }
'"$(cat "$system.potential")"'
!/^#/ {
  r = $1 ; e = pot(r) ; wgt = r*r * exp(-e/(gasR*temp)) ;
  n = n+1 ; rv[n] = r ; f[n] = wgt ; swgt += wgt ;
  }
END {
  rs = 0 ; ps = 0 ;
  scl = swgt * (rv[n]-rv[1])/(n-1) ;
  for (i=1;i<=n;i++) {
     printf " %10.5f %20.15f\n", rv[i], f[i]/scl ;
     if (i<n) {
       ps += (rv[i+1]-rv[i])*(f[i]/scl) ;
       rs += (rv[i+1]+rv[i])*(rv[i+1]-rv[i])*(f[i]/scl)/2 ;
       }
     }
  printf "Classical: total probability = %12.6f, <r> = %12.6f\n", ps, rs > "/dev/stderr" ;
  }
' temp="$temp" > "$system.classical.$temp.dat"
#
#  Plot radial distribution
#
( tee "$system.plot.$temp.gpl" | gnuplot -persist ) <<EOI
set grid
plot "$system.classical.$temp.dat" title "classical, T=$temp" with linespoints 1, \
     "$system.quantum.$temp.dat"   title "quantum, T=$temp"   with linespoints 2
EOI
#
#  Plot levels populations:
#
awk < $system.quantum.$temp.dat > $system.populations.$temp.dat '
BEGIN { maxl = -1 ; }
/^# /&&(NF==4) {
  l = $2 + 0 ; n = $3 + 0 ; pop = $4 + 0 ;
  if (l>maxl) { maxl = l ; }
  totpop[l] += $4 ;
  }
END {
  printf "# L  Population\n" ;
  for (l=0;l<=maxl;l++) {
    printf " %3d %10.7f\n", l, totpop[l] ;
    }
  }' ;
( tee $system.populations.$temp.gpl | gnuplot -persist ) << EOI
set grid
set xlabel "L"
set ylabel "Population"
plot "$system.populations.$temp.dat" title "T=$temp" with linespoints
EOI
