cmake_minimum_required( VERSION 2.6 )

PROJECT( BIMAP )

find_package(MPI)

add_subdirectory( src )

enable_testing()

add_subdirectory(test)
