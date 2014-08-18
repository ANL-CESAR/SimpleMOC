#!/bin/bash
for i in {1..32}
do
	./SimpleMOC -i default.in -t $i
done
