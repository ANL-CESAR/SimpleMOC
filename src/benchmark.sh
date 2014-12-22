#!/bin/bash
for i in {1..32}
do
	./SimpleMOC -i big.in -t $i
done
