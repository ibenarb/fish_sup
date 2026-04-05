#!/usr/bin/env bash
set -euo pipefail

cmake -S . -B build
cmake --build build
./build/fish_sup 2 3 30 40 30 20 20 10 false
