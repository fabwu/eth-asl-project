#!/bin/bash

REPETITIONS="1"
dir=$(dirname "$(readlink -f "$0")")
cd "$dir" || exit 1

for img in lena_64 monkey_128 lena_256 lena_512 grey-parrot_1024 matterhorn_2048
do
  echo " $img"
  echo "--------------------------------"
  # for exe in 09_basel 08_remove_function_and_add_fma 05_rtd_no_precompute
  for exe in 01_precompute_indices
  do
    echo "[$img - $exe] start (`date +"%Y-%m-%dT %H:%M:%S"`)"
    mkdir -p "${dir}/data/${exe}"
    ${dir}/../../build/${exe} \
                  -b -c -i 10 \
                  -r "$REPETITIONS" \
                  -f "${dir}/../../images/${img}.gray" \
                  -o "${dir}/data/${exe}/${img}.csv"
    echo "--------------------------------"
  done
done
