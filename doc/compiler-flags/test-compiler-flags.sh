#!/bin/bash

dir=$(dirname "$(readlink -f "$0")")
cd "$dir/../../" || exit 1

echo "----------------------------"
echo "  set no flags"
echo "----------------------------"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DOPT_ARCH=NO \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=NO \
      -DOPT_OFAST=NO \
      -DOPT_VEC=NO \
      -DOPT_UNSAFE_MATH=NO \
      -DOPT_UNROLL_LOOPS=NO \
      -DOPT_FAST_MATH=NO \
      -B build/
cmake --build build --target 37_full_simd
./build/37_full_simd -b -c -r 1 -f images/lena_512.gray


echo "----------------------------"
echo "  add O3"
echo "----------------------------"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DOPT_ARCH=NO \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=YES \
      -DOPT_OFAST=NO \
      -DOPT_VEC=NO \
      -DOPT_UNSAFE_MATH=NO \
      -DOPT_UNROLL_LOOPS=NO \
      -DOPT_FAST_MATH=NO \
      -B build/
cmake --build build --target 37_full_simd
./build/37_full_simd -b -c -r 1 -f images/lena_512.gray

# fast but  still stable
echo "----------------------------"
echo "  now also vectorize"
echo "----------------------------"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DOPT_ARCH=YES \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=YES \
      -DOPT_OFAST=NO \
      -DOPT_VEC=YES \
      -DOPT_UNSAFE_MATH=NO \
      -DOPT_UNROLL_LOOPS=NO \
      -DOPT_FAST_MATH=NO \
      -B build/
cmake --build build --target 37_full_simd
./build/37_full_simd -b -c -r 1 -f images/lena_512.gray

echo "----------------------------"
echo "  now also unroll loop"
echo "----------------------------"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DOPT_ARCH=YES \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=YES \
      -DOPT_OFAST=NO \
      -DOPT_VEC=YES \
      -DOPT_UNSAFE_MATH=NO \
      -DOPT_UNROLL_LOOPS=YES \
      -DOPT_FAST_MATH=YES \
      -B build/

echo "----------------------------"
echo "  set fast flags"
echo "----------------------------"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DOPT_ARCH=YES \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=YES \
      -DOPT_OFAST=NO \
      -DOPT_VEC=YES \
      -DOPT_UNSAFE_MATH=YES \
      -DOPT_UNROLL_LOOPS=YES \
      -DOPT_FAST_MATH=YES \
      -B build/

cmake --build build --target 37_full_simd
./build/37_full_simd -b -c -r 1 -f images/lena_512.gray


echo "----------------------------"
echo "  set fast flags"
echo "----------------------------"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DOPT_ARCH=YES \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=YES \
      -DOPT_OFAST=YES \
      -DOPT_VEC=YES \
      -DOPT_UNSAFE_MATH=YES \
      -DOPT_UNROLL_LOOPS=YES \
      -DOPT_FAST_MATH=YES \
      -B build/

cmake --build build --target 37_full_simd
./build/37_full_simd -b -c -r 1 -f images/lena_512.gray
