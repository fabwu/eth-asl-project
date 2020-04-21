# ASL - TEAM 002

Usage example:

```
# configure build type and specify optimazations
cmake -DCMAKE_BUILD_TYPE=Release -DOPT_ARCH=YES -B build/

# build e.g. benchmark binary
cmake --build build --target benchmark_compress

# run benchmark
./build/benchmark_compress images/lenasmall.gray
```
