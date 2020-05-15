# ASL - TEAM 002

Usage example:

```
# configure build type and specify optimazations
cmake -DCMAKE_BUILD_TYPE=Release -DOPT_ARCH=YES -B build/

# build e.g. benchmark binary
cmake --build build --target fic

# run benchmark
./build/fic -b -c -f images/lenasmall.gray
```

# Git LFS

Files in the images folder are stored with git lfs.

Install git lfs on your host and then run `git lfs install` and `git lfs pull`
here.
