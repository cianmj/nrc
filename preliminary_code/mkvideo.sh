#!/bin/sh
./gnup.sh $1
cd images
./video-convert.sh
./video-compress.sh $2
exit 0
