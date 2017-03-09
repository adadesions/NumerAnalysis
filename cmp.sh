#/bin/bash
FILENAME=$1
echo "=== Start Processing ==="
g++ $FILENAME -o $FILENAME.o 
./$FILENAME.o
echo "=== Done Processing ==="
echo "By Ada"
