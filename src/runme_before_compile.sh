#!/bin/bash
MASTER=thunder-master

for i in include Makefile MACHINES a.GLOBAL b.FUNCTIONS c.SYSTEM g.XC_FUNCTIONALS p.THEORY 
do
    if [ -e $i ]
    then
	rm -rf $i
    fi
    ln -s ../../${MASTER}/src/$i
done
