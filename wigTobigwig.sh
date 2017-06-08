#!/bin/sh
/home/tempGC/programs/wigToBigWig $1 $2 $(dirname $1)/$(basename $1 .wig.gz).bw