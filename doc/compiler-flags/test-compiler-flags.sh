#!/bin/bash

REPETITIONS=2
images=("lion_64" "lion_128" "lion_256" "lion_512" "lion_1024" "lion_2048" "lion_4096")

# test images
# images=("lion_64")

dir=$(dirname "$(readlink -f "$0")")
cd "$dir/../../" || exit 1
rm -rf ${dir}/../../build

runbenchmarks40(){
      cmake --build build --target 25_ilp
      name="$1"
      exe="25_ilp"
      for img in ${images[@]}
      do
            echo "--------------------------------"
            echo "[$img - $exe] start (`date +"%Y-%m-%dT %H:%M:%S"`)"
            mkdir -p "${dir}/data/${name}-${exe}"
            ${dir}/../../build/${exe} \
                  -b -c -i 10 \
                  -r "$REPETITIONS" \
                  -f "${dir}/../../images/${img}.gray" \
                  -o "${dir}/data/${name}-${exe}/${img}.csv"
      done
}


runbenchmarks41(){
      cmake --build build --target 51_simd_improved_rot
      name="$1"
      exe="51_simd_improved_rot"
      for img in ${images[@]}
      do
            echo "--------------------------------"
            echo "[$img - $exe] start (`date +"%Y-%m-%dT %H:%M:%S"`)"
            mkdir -p "${dir}/data/${name}-${exe}"
            ${dir}/../../build/${exe} \
                  -b -c -i 10 \
                  -r "$REPETITIONS" \
                  -f "${dir}/../../images/${img}.gray" \
                  -o "${dir}/data/${name}-${exe}/${img}.csv"
      done
}

echo "----------------------------"
echo "  gcc_Ofast"
echo "----------------------------"
rm -rf ${dir}/../../build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DICC=NO \
      -DOPT_ARCH=YES \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=NO \
      -DOPT_OFAST=YES \
      -DOPT_FMA=YES \
      -DOPT_UNROLL_LOOPS=NO \
      -B build/
runbenchmarks40 "gcc_Ofast"
# runbenchmarks41 "gcc_Ofast"


echo "----------------------------"
echo "  gcc_03_unroll "
echo "----------------------------"
rm -rf ${dir}/../../build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DICC=NO \
      -DOPT_ARCH=YES \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=YES \
      -DOPT_OFAST=NO \
      -DOPT_FMA=YES \
      -DOPT_UNROLL_LOOPS=YES \
      -B build/
runbenchmarks40 "gcc_O3_unroll"
# runbenchmarks41 "gcc_O3_unroll"


echo "----------------------------"
echo "  gcc_O3"
echo "----------------------------"
rm -rf ${dir}/../../build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DICC=NO \
      -DOPT_ARCH=YES \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=YES \
      -DOPT_OFAST=NO \
      -DOPT_FMA=YES \
      -DOPT_UNROLL_LOOPS=NO \
      -B build/
runbenchmarks40 "gcc_O3"
# runbenchmarks41 "gcc_O3"


echo "----------------------------"
echo "  gcc_03_no_fma"
echo "----------------------------"
rm -rf ${dir}/../../build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DICC=NO \
      -DOPT_ARCH=YES \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=YES \
      -DOPT_OFAST=NO \
      -DOPT_FMA=NO \
      -DOPT_UNROLL_LOOPS=NO \
      -B build/
runbenchmarks40 "gcc_O3_no_fma"
# no fma runbenchmarks41 "gcc_O3"


# echo "############################"
# echo "############################"
# echo "         ICC"
# echo "############################"
# echo "############################"


echo "----------------------------"
echo "  icc_Ofast"
echo "----------------------------"
rm -rf ${dir}/../../build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DICC=YES \
      -DOPT_ARCH=YES \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=NO \
      -DOPT_OFAST=YES \
      -DOPT_FMA=YES \
      -DOPT_UNROLL_LOOPS=NO \
      -B build/
runbenchmarks40 "icc_Ofast"
# runbenchmarks41 "icc_Ofast"


echo "----------------------------"
echo "  icc_03"
echo "----------------------------"
rm -rf ${dir}/../../build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DICC=YES \
      -DOPT_ARCH=YES \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=YES \
      -DOPT_OFAST=NO \
      -DOPT_FMA=YES \
      -DOPT_UNROLL_LOOPS=NO \
      -B build/
runbenchmarks40 "icc_O3"
# runbenchmarks41 "icc_O3"


###############################################
# Slower
###############################################

echo "----------------------------"
echo "  gcc_O2"
echo "----------------------------"
rm -rf ${dir}/../../build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DICC=NO \
      -DOPT_ARCH=YES \
      -DOPT_O1=NO \
      -DOPT_O2=YES \
      -DOPT_O3=NO \
      -DOPT_OFAST=NO \
      -DOPT_FMA=YES \
      -DOPT_UNROLL_LOOPS=NO \
      -B build/
runbenchmarks40 "gcc_O2"
# runbenchmarks41 "gcc_O2"

echo "----------------------------"
echo "  gcc_O1"
echo "----------------------------"
rm -rf ${dir}/../../build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DICC=NO \
      -DOPT_ARCH=YES \
      -DOPT_O1=YES \
      -DOPT_O2=NO \
      -DOPT_O3=NO \
      -DOPT_OFAST=NO \
      -DOPT_FMA=YES \
      -DOPT_UNROLL_LOOPS=NO \
      -B build/
runbenchmarks40 "gcc_O1"
# runbenchmarks41 "gcc_O1"

echo "----------------------------"
echo "  gcc"
echo "----------------------------"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DICC=NO \
      -DOPT_ARCH=YES \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=NO \
      -DOPT_OFAST=NO \
      -DOPT_FMA=YES \
      -DOPT_UNROLL_LOOPS=NO \
      -B build/
runbenchmarks40 "gcc"
# runbenchmarks41 "gcc"
