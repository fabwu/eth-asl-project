#!/bin/bash

dir=$(dirname "$(readlink -f "$0")")
cd "$dir"
echo "$dir"

for img in lena_64 monkey_128 lena_256 lena_512 grey-parrot_1024
do
  echo "$img start (`date +"%Y-%m-%dT %H:%M:%S"`)"
  ./../../build/fic -b -c -i 10 \
                    -f "${dir}/../../images/${img}.gray" \
                    -o "${dir}/data/${img}.csv"
  echo "$img end (`date +"%Y-%m-%dT %H:%M:%S"`)"
  echo "--------------------------------"
done
