#!/bin/bash

MASTER=fireball-master

for i in a.GLOBAL b.FUNCTIONS c.SYSTEM f.MPI g.XC_FUNCTIONALS j.ASSEMBLERS p.THEORY o.OUTPUT include Makefile.in MACHINES libs
do
    if [ -e $i ]
    then
	rm -rf $i
    fi
    ln -s ../../${MASTER}/src/$i
done
