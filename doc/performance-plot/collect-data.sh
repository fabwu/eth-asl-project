#!/bin/bash

dir=$(dirname "$(readlink -f "$0")")
cd "$dir" || exit 1


REPETITIONS="3"
for img in lion_64 lion_128 lion_256 lion_512 lion_1024
do
  echo " $img"
  echo "--------------------------------"
  for exe in 41_simd_norot_90_270 40_ilp_norot_90_270 35_simd 31_simd_precomp_rotations_no_bac_simd 25_ilp
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
for img in lion_64 lion_2048 lion_4096
do
  echo " $img"
  echo "--------------------------------"
  for exe in 41_simd_norot_90_270 40_ilp_norot_90_270 35_simd 31_simd_precomp_rotations_no_bac_simd 25_ilp
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
for img in lion_64 lion_128 lion_256 lion_512 lion_1024 lion_2048 lion_4096
do
  echo " $img"
  echo "--------------------------------"
  for exe in 0_baseline
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
