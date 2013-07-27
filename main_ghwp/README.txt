 ! 
 !  Time propagation of a 1D Gauss-Hermite wavepacket
 !
 !  Last modified: August 17, 2004
 !
 !  Written by: S.P and C.M.-J.
 !
 !  The basis function form is:
 !
 !                      1          / 2 A \1/4   /     1/2       \
 !    \phi_n(x) = ---------------- | --- |    H | (2A)   (x-R) )|  
 !                (2**n n!)**(1/2) \ Pi  /     n\               /
 !
 !                   /           \    /              2 \
 !                exp| i P*(x-R) | exp| (-A+iB) (x-R)  |
 !                   \           /    \                /
 !
 
 ghwp1D.f90
 
 
