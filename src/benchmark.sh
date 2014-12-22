#!/bin/bash
for i in 1 2 4 8 16 24 32
do
	./SimpleMOC -i big.in -t $i
done
