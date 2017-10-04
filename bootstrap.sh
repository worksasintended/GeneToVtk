#!/bin/bash

rm -r CMakeFiles
rm CMakeCache.txt

CC=gcc CXX=g++ cmake $FLAGS .
