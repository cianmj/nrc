#!/bin/sh
./gnup.sh $1
cd videos 
./video-convert.sh
./video-compress.sh $2
exit 0
