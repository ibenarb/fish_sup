cmake_minimum_required(VERSION 3.28)
project(fish_sup LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(fish_sup
        src/main.cpp
        src/types.h
        src/two_by_c.h
        src/combinatorics.h
        src/real.h
        src/enumeration.h
        src/table_record.h
        src/table_collection.h
        src/minlike.h
        src/probability.h
        src/event.h
        src/event_logic.h
        src/event_collection.h
        src/logsumexp.h
)

target_compile_options(fish_sup PRIVATE
        -Wall
        -Wextra
        -Wpedantic
)

target_compile_options(fish_sup PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-DNDEBUG>
)

target_link_libraries(fish_sup PRIVATE
        gmp
        gmpxx
)

target_link_options(fish_sup PRIVATE
        $<$<CONFIG:Release>:-s>
)