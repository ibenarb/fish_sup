# Build on Linux / WSL / Ubuntu

## Requirements

- CMake >= 3.28
- C++20 compiler
- GMP development package

## Install dependencies on Ubuntu

```text
sudo apt update && sudo apt install -y build-essential cmake libgmp-dev
```

## Configure and build

```text
cmake -S . -B build && cmake --build build
```

## Run the reference example

```text
./build/fish_sup 2 3 30 40 30 20 20 10 false
```
