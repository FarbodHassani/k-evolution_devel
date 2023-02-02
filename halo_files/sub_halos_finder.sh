#!/bin/bash

# Set boxsize variable
boxsize=1200.000 # In the rockstar makefile you need to set DRADIUS_CONVERSION=1.0 for the parents: part

find . -name "*.list" -exec sh -c './find_parents "$0" boxsize > "${0%.*}_sub-halo.list" ' {} \;
