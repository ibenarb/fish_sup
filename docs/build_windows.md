# Build on Windows

## Requirements

- CMake (>= 3.28)
- A C++20 capable compiler (e.g. MinGW g++)
- GMP and gmpxx libraries

## Notes on setup

The exact setup depends on your Windows toolchain.

Typical options include:

- MinGW-w64 with g++
- MSYS2 environments
- CLion with bundled MinGW toolchain

Ensure that GMP and gmpxx are available to the compiler and linker.

## Configure and build (MinGW / MSYS2 example)

```text
cmake -S . -B build
cmake --build build
```

## Run the program

```text
build\\fish_sup.exe 2 3 30 40 30 20 20 10 false
```

## Remarks

Floating-point output should match the Linux reference up to final rounding differences if the same precision settings are used.
