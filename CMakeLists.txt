cmake_minimum_required(VERSION 3.20)
project(DNAMotifFinder VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++23")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wconversion -Wshadow")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter -Wno-sign-conversion")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmodules-ts -fconcepts-diagnostics-depth=3")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wconversion -Wshadow")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter -Wno-sign-conversion")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmodules -fconcepts-diagnostics-depth=3")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /permissive- /std:c++23")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /experimental:module")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -DDEBUG -fsanitize=address,undefined")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -march=native -flto")

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(include)
include_directories(src)

set(SOURCES
    src/main.cpp
    src/dna_parser.cpp
    src/motif_finder.cpp
    src/iupac_codes.cpp
    src/mpi_manager.cpp
    src/parallel_processor.cpp
)

set(HEADERS
    include/dna_parser.h
    include/motif_finder.h
    include/iupac_codes.h
    include/mpi_manager.h
    include/parallel_processor.h
    include/common.h
    include/concepts.h
)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 23
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)

target_link_libraries(${PROJECT_NAME}
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
)

find_package(GTest REQUIRED)

enable_testing()

add_subdirectory(tests)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
install(FILES ${HEADERS} DESTINATION include)

add_custom_target(format
    COMMAND find ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/include -name "*.cpp" -o -name "*.h" | xargs clang-format -i
    COMMENT "Formatting source code"
)

add_custom_target(analyze
    COMMAND cppcheck --enable=all --std=c++23 --suppress=missingIncludeSystem ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/include
    COMMENT "Running static analysis"
)
