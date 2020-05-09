# Valgrind

It is important to compile with debug options.

## Cachegrind

```
valgrind --tool=cachegrind ./build/fic -c -d -e 100 -f images/lena_256.gray
```

## Callgrind

```
valgrind --tool=callgrind ./build/fic -c -d -e 100 -f images/lena_256.gray --log-file="callgrind
```
