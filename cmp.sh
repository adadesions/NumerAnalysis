#!/bin/bash
FILENAME=$1
g++ $FILENAME -o $FILENAME.o 
./$FILENAME.o
echo "=== Done Processing ==="
