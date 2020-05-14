#!/bin/bash


dir=$(dirname "$(readlink -f "$0")")
cd "$dir" || exit 1

for exe in baseline 01_precompute_indices 02_precomputations 03_jonas
do
  echo " $exe"
  echo "--------------------------------"
  for img in lena_64 monkey_128 lena_256 lena_512 grey-parrot_1024
  do
    echo "$img start (`date +"%Y-%m-%dT %H:%M:%S"`)"
    mkdir -p "${dir}/data/${exe}"
    ${dir}/../../build/${exe} \
                  -b -c -i 10 \
                  -f "${dir}/../../images/${img}.gray" \
                  -o "${dir}/data/${exe}/${img}.csv"
    echo "$img end (`date +"%Y-%m-%dT %H:%M:%S"`)"
    echo "--------------------------------"
  done
done
