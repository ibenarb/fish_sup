#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
SRC_FILE="$PROJECT_ROOT/src/main.cpp"
LINUX_OUT="$PROJECT_ROOT/fish_sup"
WIN_OUT="$PROJECT_ROOT/fish_sup.exe"
MINGW_PREFIX="/mnt/c/msys64/mingw64"
MINGW_INCLUDE="$MINGW_PREFIX/include"
MINGW_LIB_GMPXX="$MINGW_PREFIX/lib/libgmpxx.a"
MINGW_LIB_GMP="$MINGW_PREFIX/lib/libgmp.a"
LINUX_CXX="${CXX:-g++}"
WIN_CXX="x86_64-w64-mingw32-g++"
COMMON_FLAGS=(-std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic)
echo "==> Building Linux binary"
"$LINUX_CXX" "${COMMON_FLAGS[@]}" "$SRC_FILE" -lgmpxx -lgmp -o "$LINUX_OUT"
echo "==> Building Windows binary"
"$WIN_CXX" "${COMMON_FLAGS[@]}" -I"$MINGW_INCLUDE" "$SRC_FILE" "$MINGW_LIB_GMPXX" "$MINGW_LIB_GMP" -o "$WIN_OUT" -static-libstdc++ -static-libgcc
echo
echo "Build finished:"
echo "  Linux   -> $LINUX_OUT"
echo "  Windows -> $WIN_OUT"
