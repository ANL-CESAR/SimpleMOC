#!/bin/bash
for i in {1..32}
do
	./SimpleMOC -i bgq.in -t $i
done
