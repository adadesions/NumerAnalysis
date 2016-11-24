#!/bin/bash
FILENAME=$1
g++ $FILENAME.cpp -o $FILENAME.out
./$FILENAME.out
