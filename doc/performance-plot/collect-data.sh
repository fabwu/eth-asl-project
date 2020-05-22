#!/bin/bash

dir=$(dirname "$(readlink -f "$0")")
cd "$dir" || exit 1

REPETITIONS="3"
for img in lena_64 lena_512 grey-parrot_1024
do
  echo " $img"
  echo "--------------------------------"
  for exe in 36_pascal 35_simd 31_simd_precomp_rotations_no_bac_simd 30_simd_precomp_rotations 25_ilp
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

REPETITIONS="1"
for img in monkey_128 matterhorn_2048 lion_4096
do
  echo " $img"
  echo "--------------------------------"
  for exe in 36_pascal 35_simd 31_simd_precomp_rotations_no_bac_simd 30_simd_precomp_rotations 25_ilp
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
