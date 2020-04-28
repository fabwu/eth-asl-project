#!/bin/bash

./build/fic -b -c -i 10 -s 8 16 -f images/lena_64.gray \
            -o doc/performance-plot/data/lena_64.csv
./build/fic -b -c -i 10 -s 8 16 -f images/monkey_128.gray \
            -o doc/performance-plot/data/monkey_128.csv
./build/fic -b -c -i 10 -s 8 16 -f images/lena_256.gray \
            -o doc/performance-plot/data/lena_256.csv
./build/fic -b -c -i 10 -s 8 16 -f images/lena_512.gray \
            -o doc/performance-plot/data/lena_512.csv
./build/fic -b -c -i 10 -s 8 16 -f images/grey-parrot_1024.gray \
            -o doc/performance-plot/data/grey-parrot_1024.csv
