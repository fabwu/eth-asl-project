#!/bin/bash

REPETITIONS=2
images=("lion_64" "lion_128" "lion_256" "lion_512" "lion_1024" "lion_2048" "lion_4096")

# test images
# images=("lion_64" "lion_128")

dir=$(dirname "$(readlink -f "$0")")
cd "$dir/../../" || exit 1

runbenchmarks40(){
      cmake --build build --target 40_ilp_norot_90_270
      name="$1"
      exe="40_ilp_norot_90_270"
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
      cmake --build build --target 41_simd_norot_90_270
      name="$1"
      exe="41_simd_norot_90_270"
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
echo "  gcc_Ofast_fma"
echo "----------------------------"
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
runbenchmarks40 "gcc_Ofast_fma"
runbenchmarks41 "gcc_Ofast_fma"


echo "----------------------------"
echo "  gcc_03_fma_unroll "
echo "----------------------------"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DICC=NO \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=YES \
      -DOPT_OFAST=NO \
      -DOPT_FMA=YES \
      -DOPT_UNROLL_LOOPS=YES \
      -B build/
runbenchmarks40 "gcc_O3_fma_unroll"
runbenchmarks41 "gcc_O3_fma_unroll"


echo "----------------------------"
echo "  gcc_O3_fma"
echo "----------------------------"
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
runbenchmarks40 "gcc_O3_fma"
runbenchmarks41 "gcc_O3_fma"


echo "----------------------------"
echo "  gcc_03"
echo "----------------------------"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DICC=NO \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=YES \
      -DOPT_OFAST=NO \
      -DOPT_FMA=NO \
      -DOPT_UNROLL_LOOPS=NO \
      -B build/
runbenchmarks40 "gcc_O3"
# no fma runbenchmarks41 "gcc_O3"


echo "############################"
echo "############################"
echo "         ICC"
echo "############################"
echo "############################"


echo "----------------------------"
echo "  icc_Ofast_fma"
echo "----------------------------"
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
runbenchmarks41 "icc_Ofast_fma"
runbenchmarks40 "icc_Ofast_fma"


echo "----------------------------"
echo "  icc_03"
echo "----------------------------"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DICC=NO \
      -DOPT_O1=NO \
      -DOPT_O2=NO \
      -DOPT_O3=YES \
      -DOPT_OFAST=NO \
      -DOPT_FMA=NO \
      -DOPT_UNROLL_LOOPS=NO \
      -B build/
runbenchmarks40 "icc_O3"
# no fma runbenchmarks41 "icc_O3"
