#!/bin/sh
echo "search for... "
read word
echo " ... "
grep -i -n "$word" *.f90
echo "done search"
exit 0
