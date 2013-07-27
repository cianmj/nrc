#!/bin/sh
gprof run.x gmon.out > profile.out
rm gmon.out
exit 0
