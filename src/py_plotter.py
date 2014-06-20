#!/usr/bin/env python

import matplotlib.pyplot as plt
import random

with open("gen_points_2D.txt") as fh:
	lines = fh.readlines()

# store data in vectors
x = []
y = []
z = []
assembly = []
pin = []
zone = []
index = []

# option to only read a fraction of the data
fraction = 1

# skipping header line, reading the rest
for i in range( len(lines)/fraction - 1 ):
	words = lines[i+1].split()
	x.append( float(words[0]) )
	y.append( float(words[1]) )
	z.append( float(words[2]) )
	assembly.append ( int(words[3]) )
	pin.append( int(words[4]) )
	zone.append( int(words[5]) )
	index.append( int(words[6]) )

print "Finished loading data..."

# define which attribute to color by
attribute = pin

# calculate appropriate marker size
markerSize = 5 * 10**6 / float(len(attribute))

# create array of random numbers with a length of the maximum attribute
max_attribute = max(attribute)
randoms = []
for i in range(max_attribute + 1):
	randoms.append(random.random())
	
colors = []
for i in range( len(attribute) ):
	colors.append( randoms[attribute[i]] )

print "Finished calculating colors..."

print "Marker Size = " + str(markerSize)

plt.scatter(x, y, marker='.', s=markerSize, linewidths=0, c=colors, cmap=plt.cm.Set1)

plt.show()
