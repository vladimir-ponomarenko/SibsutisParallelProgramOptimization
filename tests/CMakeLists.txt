cmake_minimum_required(VERSION 3.16)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GTest REQUIRED)

set(TEST_SOURCES
    test_iupac_codes.cpp
    test_dna_parser.cpp
    test_motif_finder.cpp
    test_mpi_simple.cpp
    test_main.cpp
    ../src/iupac_codes.cpp
    ../src/dna_parser.cpp
    ../src/motif_finder.cpp
    ../src/mpi_manager.cpp
)

add_executable(dna_motif_tests ${TEST_SOURCES})

target_include_directories(dna_motif_tests PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(dna_motif_tests 
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
    GTest::GTest
    GTest::Main
)

add_test(NAME iupac_codes_test COMMAND dna_motif_tests --gtest_filter=IUPACCodesTest.*)
add_test(NAME dna_parser_test COMMAND dna_motif_tests --gtest_filter=DNAParserTest.*)
add_test(NAME motif_finder_test COMMAND dna_motif_tests --gtest_filter=MotifFinderTest.*)
add_test(NAME mpi_simple_test COMMAND dna_motif_tests --gtest_filter=MPIManagerSimpleTest.*)
add_test(NAME main_test COMMAND dna_motif_tests --gtest_filter=MainTest.*)

set_tests_properties(iupac_codes_test PROPERTIES TIMEOUT 30)
set_tests_properties(dna_parser_test PROPERTIES TIMEOUT 30)
set_tests_properties(motif_finder_test PROPERTIES TIMEOUT 30)
set_tests_properties(mpi_simple_test PROPERTIES TIMEOUT 30)
set_tests_properties(main_test PROPERTIES TIMEOUT 30)
