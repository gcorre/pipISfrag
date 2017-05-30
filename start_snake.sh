#!/bin/sh
snakemake -s /home/tempGC/Scripts/snakefile2.hs -d $PWD --dag | dot -Tpdf > dag.pdf;
snakemake -s /home/tempGC/Scripts/snakefile2.hs -j 8 -d $PWD > pipeline.log 2>pip.err;
snakemake -s /home/tempGC/Scripts/snakefile2.hs -D -d $PWD > files.log;