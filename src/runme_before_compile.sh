#!/bin/bash
MASTER=thunder-master

for i in include a.GLOBAL b.FUNCTIONS c.SYSTEM g.XC_FUNCTIONALS p.THEORY Makefile MACHINES
do
    if [ -e $i ]
    then
	rm -rf $i
    fi
    ln -s ../../${MASTER}/src/$i
done
