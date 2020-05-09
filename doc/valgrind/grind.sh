#!/bin/bash

dir=$(dirname "$(readlink -f "$0")")
cd "$dir/../../" || exit 1
mkdir -p "$dir/logs"

fic_opts="-c -d -e 100 -f images/lena_256.gray"

# baseline
valgrind --tool=cachegrind \
         --cachegrind-out-file="./doc/valgrind/logs/cachegrind-baseline.log" \
         ./build/baseline $fic_opts &
valgrind --tool=callgrind \
         --callgrind-out-file="./doc/valgrind/logs/callgrind-baseline.log" \
         ./build/baseline $fic_opts &

# 01_precompute_indices
valgrind --tool=cachegrind \
         --cachegrind-out-file="./doc/valgrind/logs/cachegrind-01_precompute_indices.log" \
         ./build/01_precompute_indices $fic_opts &
valgrind --tool=callgrind \
         --callgrind-out-file="./doc/valgrind/logs/callgrind-01_precompute_indices.log" \
         ./build/01_precompute_indices $fic_opts &
